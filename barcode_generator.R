#!/usr/bin/env Rscript

library("gtools")
library("tidyverse")

source("config.R")

# initializing empty variables

bc.length <- integer(iter.number)
optimal.barcodes <- NULL
optimal.dist <- NULL

#############################
## Initial Barcode Filtering
#############################

## Generate Barcodes
## ----------

bases <- c("A", "T", "C", "G")
barcodes <- permutations(4, barcode.length, bases, repeats.allowed=T) %>%
	as_tibble %>%
	unite(col="barcodes", sep="")

print(paste(nrow(barcodes), "initial barcodes"))

## Remove Truseq Adapters
## ----------

## load truseq adapters

truseq <- read.delim("truseq_barcodes.tsv", header=F, sep="\t", stringsAsFactors=F) %>%
	as_tibble %>%
	dplyr::rename("barcodes"=1)

## remove truseq adapters from barcodes

barcodes <- truseq %>%
	mutate("truseq"="Y") %>%
	right_join(., barcodes, by="barcodes")

remaining <- filter(barcodes, is.na(truseq))
print(paste(nrow(remaining), "barcodes after removal of truseq barcodes"))

## Filter Barcodes by GC Content
## ----------

## keep barcodes with GC content >= 0.25 and <= 0.75

barcodes <- mutate(barcodes, "GC"= str_count(barcodes, pattern=c("G","C")) / barcode.length)

remaining <- filter(barcodes, GC <= 0.75 & GC >= 0.25)
print(paste(nrow(remaining), "barcodes after GC content filter"))

## Filter Identical Base Runs
## ----------

## remove barcodes with 4+ runs of identical bases

barcodes <- mutate(
	barcodes,
	single_base_runs = case_when(grepl(barcodes, pattern="(A{4,}|T{4,}|G{4,}|C{4,})") ~ "Y")
)

remaining <- filter(barcodes, is.na(single_base_runs))
print(paste(nrow(remaining), "barcodes after 4+ sequential base filter"))

## Remove Barcodes Similar to Truseq
## ----------

## remove barcodes that are less than n.dist away from a truseq barcode

barcodes <- adist(pull(barcodes, barcodes), pull(truseq, barcodes)) %>%
	as_tibble %>%
	mutate_all(~ . < n.dist) %>%
	rowSums %>%
	map_chr(~ifelse(.==1, "Y", NA)) %>%
	as_tibble %>%
	rename("truseq_similarity"=1) %>%
	mutate("barcodes"=pull(barcodes, barcodes)) %>%
	left_join(barcodes, ., by="barcodes")

remaining <- filter(remaining, is.na(truseq_similarity))
print(paste(length(remaining), "barcodes after removing barcodes with a dist <", n.dist, "from truseq barcodes"))

########################################
## Maximize Barcode Number and Distance
########################################

for (iter in 1:iter.number) {
print(paste(".....iteration number:", paste0(iter,"/",iter.number), "....."))

## jumble up barcode order and remove barcodes that are too similar (distance < n.dist)
## ----------

## jumble barcode order

barcode.sample <- remaining %>%
	pull(barcodes) %>%
	sample

## get barcode distances

barcode.dists <- adist(barcode.sample, barcode.sample)
rownames(barcode.dists) <- barcode.sample
colnames(barcode.dists) <- barcode.sample

## keep only one sequence from sequences with distances < n.dist

discarded.barcodes <- c()
for (n in 1:nrow(barcode.dists)) {
	# skip if barcode was already discarded
	if (colnames(barcode.dists)[n] %in% discarded.barcodes) next

	# get barcodes with too small of a distance
	to.discard <- which(barcode.dists[n,] < n.dist & barcode.dists[n,] > 0) %>%
		colnames(barcode.dists)[.] %>%
		discard(. %in% discarded.barcodes)

	# add barcodes with too small of a distance to exclusion list
	discarded.barcodes <- c(discarded.barcodes, to.discard)
}

barcode.sample <- discard(barcode.sample, barcode.sample %in% discarded.barcodes)

print(paste(length(barcode.sample), "barcodes after keeping only one barcode of those with distances <", n.dist))

## remove reverse complements and palindromes
## ----------

## function to reverse complement a sequence

reverse.complement <- function(x) {
	x %>%
		str_split(pattern="", simplify=TRUE) %>%
		rev %>%
		map_chr(
			~ if (. == "A") {
				return("T")
			} else if (. == "T") {
				return("A")
			} else if (. == "G") {
				return("C")
			} else {
				return("G")
			}
		) %>%
		paste(., collapse="") %>%
		return
}

## check the reverse complement barcodes for any matches in the barcode list

removed.complements <- c()
temp <- map_chr(
	barcode.sample,
	~ if (. %in% removed.complements) {
		return(NA)
	} else {
		rev.comp <- reverse.complement(.)
		if (rev.comp %in% barcode.sample) {
			removed.complements <<- c(removed.complements, rev.comp)
			return(rev.comp)
		} else {
			return(NA)
		}
	}
)

## remove the barcodes that had reverse complements and palindromes

barcode.sample <- discard(barcode.sample, barcode.sample %in% removed.complements)

print(paste(length(barcode.sample), "barcodes after removing barcodes that are reverse complements or palindromes"))


## Reduce Position Matching 3mers
## ----------

# using a window size of 3, check if any barcodes have matching 3-mers in the same base positions
# if there are more than max.nmer matches, randomly remove matching barcodes until you are down to only max.nmer matches

keep <- barcode.sample

for (n in 1:(barcode.length-2)) {
	# pull 3mer from current position
	barcode.df <- keep %>%
		as_tibble %>%
		dplyr::rename("barcodes"=1) %>%
		mutate("slice"=substr(barcodes, n, n+2))

	# get counts for each unique kmer
	kmer.counts <- barcode.df %>%
		count(slice) %>%
		left_join(barcode.df, ., by="slice")

	# stash barcodes with 3mer already less than limit
	already.passing <- kmer.counts %>%
		filter(n <= max.nmer) %>%
		pull(barcodes)

	# for barcodes above the 3mer limit, randomly sample barcodes down to limit
	make.passing <- kmer.counts %>%
		filter(n > max.nmer) %>%
		group_by(slice) %>%
		sample_n(max.nmer, replace=FALSE) %>%
		pull(barcodes)

	# add filtered barcodes back to barcodes that were stashed earlier
	keep <- c(already.passing, make.passing)
}

kept.barcodes <- keep

print(paste(length(kept.barcodes), "after reducing the number of barcodes with more than", max.nmer, "position matching 3-mers"))


#####
# remove barcodes if a base appears more than max.base.frac of the time in one position
#####

keep <- kept.barcodes

for (n in 1:(barcode.length-2)) {
        # pull base from current position
        barcode.df <- keep %>%
                as_tibble %>%
                dplyr::rename("barcodes"=1) %>%
                mutate("base"=substr(barcodes, n, n))

        # get counts for each unique base
        base.counts <- barcode.df %>%
                count(base) %>%
		dplyr::rename("count"=n) %>%
		mutate("perc"= count / sum(count)) %>%
                left_join(barcode.df, ., by="base") 
        
	# find the number that gets you close to max.base.frac
	limits <- base.counts %>%
		count(base, perc) %>%
		dplyr::rename("count"=n) %>%
		mutate(recommended = count - floor((count - (max.base.frac * sum(count))) / (1-max.base.frac))) %>%
		left_join(base.counts, ., by=c("base", "count", "perc"))

	# store barcodes where bases aren't a problem
	already.passing <- limits %>%
		filter(count < recommended) %>%
		pull(barcodes)

	# remove random barcodes so base fractions aren't too similar
        make.passing <- limits %>%
                filter(count > recommended) %>%
		group_by(base) %>%
		group_split %>%
		map(
			~ sample_n(., pull(., recommended) %>% unique) %>%
			pull(barcodes)
		) %>%
		unlist
        
        # add filtered barcodes back to barcodes that were stashed earlier
        keep <- c(already.passing, make.passing)
}

kept.barcodes <- keep

print(paste(length(kept.barcodes), "barcodes after making sure a base doesn't appear more than", paste0(max.base.frac * 100, "%"), "in one position")) 


#####
# check if desired barcode number was found
#####


# what to do if there are enough barcodes
if (length(barcodes) >= (desired.barcodes - 24)) {
	# keep track of barcode lengths
	bc.length[iter] <- length(barcodes) + 24
 	# randomly sample desired.barcodes - 24 barcodes from barcodes list
 	print(paste(length(barcodes) - (desired.barcodes - 24), "barcodes randomly removed to achieve desired barcode number of", desired.barcodes))
 	barcodes <- sample(barcodes, desired.barcodes - 24)
 	# add truseq barcodes to list
 	barcodes <- c(truseq, barcodes)
 	print(paste(length(barcodes), "barcodes after adding back the truseq barcodes"))
 	# getting average distance between barcodes
 	dist <- adist(barcodes)
 	dist <- cbind(dist, 1:nrow(dist))
 	mean.dist <- sum(apply(dist, 1, function(x) {sum(x[-x[length(x)]])})) / ((nrow(dist)**2) - nrow(dist))
 	print(paste(mean.dist, "average distance between barcodes"))
	# check to see if the barcodes should be saved 
 	if (length(optimal.barcodes)==0) {
 		optimal.barcodes <- barcodes
 		optimal.dist <- mean.dist
	} else if (mean.dist <= optimal.dist) {
 		print(paste("the current average distance of", mean.dist, "is lower than or equal to the optimal barcode set with distance", optimal.dist))
 	} else {
 		optimal.barcodes <- barcodes
 		optimal.dist <- mean.dist
 		print(paste("the current average distance of", mean.dist, "is higher than the optimal barcode set with distance", optimal.dist))
 	}
# what to do if there are not enough barcodes
} else {
	bc.length[iter] <- length(barcodes) + 24
	print(paste(length(barcodes) + 24, "barcodes after adding back the truseq barcodes"))
	print(paste("desired barcode number of", desired.barcodes, "was not met, moving to next iteration"))
}
# what to do in last iteration
if (iter == iter.number) {
	# if number of desired barcodes was not met
	if (length(optimal.barcodes) == 0) {
		stop(paste("desired barcode number of", desired.barcodes, "could not be achieved"))
	# if number of desierd barcodes was met
	} else {
		print(".....exporting optimal barcode set.....")
  		print(optimal.barcodes)
  		#save barcodes
  		barcodes.export <- data.frame("ID"=sprintf("bc%04d", 1:desired.barcodes))
  		barcodes.export$barcodes <- optimal.barcodes
  		barcodes.export$reverse.complement <- sapply(optimal.barcodes, reverse.complement)
  		barcodes.export$adapter <- sapply(barcodes.export$reverse.complement, function(x) {paste0("CAAGCAGAAGACGGCATACGAGAT",x,"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNN")})
  		write.table(barcodes.export, "barcodes.txt", sep="\t", quote=F, col.names=T, row.names=F)
		
		bc.length.plot <- data.frame("iteration" = 1:iter.number, "bc.length" = bc.length)
		p <- ggplot(bc.length.plot, aes(bc.length)) +
			geom_histogram(fill="dodgerblue2", binwidth=1, color="white") +
			theme_classic()
  
		pdf("barcode_info.pdf")
		print(p)
		dev.off()
  		print(paste(optimal.dist, " average hamming distance for optimal barcode set"))
  		stop(paste("desired barcode number", desired.barcodes, "met, output saved"))
	}
}
}
