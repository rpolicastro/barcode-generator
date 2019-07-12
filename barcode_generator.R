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

## Filter Barcodes by GC Content
## ----------

## keep barcodes with GC content >= 0.25 and <= 0.75

barcodes <- mutate(barcodes, "GC"= str_count(barcodes, pattern=c("G","C")) / barcode.length)

discard <- barcodes %>%
	filter(GC < 0.25 | GC > 0.75) %>%
	pull(barcodes)

remaining <- barcodes %>%
	pull(barcodes) %>%
	discard(. %in% discard)
print(paste(length(remaining), "barcodes after removing barcodes with GC content <= 25% or GC content >= 75%"))

## Filter Identical Base Runs
## ----------

## remove barcodes with 4+ runs of identical bases

barcodes <- mutate(
	barcodes,
	single_base_runs = case_when(grepl(barcodes, pattern="(A{4,}|T{4,}|G{4,}|C{4,})") ~ "Y")
)

discard <- barcodes %>%
	filter(single_base_runs == "Y") %>%
	pull(barcodes)

remaining <- discard(remaining, remaining %in% discard)
print(paste(length(remaining), "barcodes after removing barcodes with 4+ identical base runs")) 

## Keep Optimal Barcodes Thus Far
## ----------

remaining <- barcodes %>%
	filter(
		GC >= 0.25 & GC <= 0.75,
		is.na(single_base_runs)
	) %>%
	pull(barcodes)

########################################
## Maximize Barcode Number and Distance
########################################

for (iter in 1:iter.number) {
print(paste(".....iteration number:", paste0(iter,"/",iter.number), "....."))

## jumble up barcode order and remove barcodes that are too similar (distance < n.dist)
## ----------

## jumble barcode order

barcode.sample <- sample(remaining, replace=FALSE)

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

## Reduce Position Matching Bases
## ----------

## remove barcodes if a base appears more than max.base.frac of the time in one position

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
		filter(count <= recommended) %>%
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

## check if desired barcode number was found
## ----------

if (length(kept.barcodes) < desired.barcodes) {
	print(paste("desired barcode number of", desired.barcodes, "was not met, moving to next iteration"))
} else {
	# sample down to desired barcode number
	kept.barcodes <- sample(kept.barcodes, desired.barcodes, replace=FALSE)
	
	# get mean distance between barcodes
	mean.dist <- kept.barcodes %>%
		adist %>%
		as.numeric %>%
		discard(. == 0) %>%
		mean
	print(paste(mean.dist, "average distance between barcodes"))

	# check to see if current barcodes should be saved
	if (length(optimal.barcodes)==0) {
		optimal.barcodes <- kept.barcodes
		optimal.dist <- mean.dist
	} else if (mean.dist <= optimal.dist) {
		print(paste("the current average distance of", mean.dist, "is lower than or equal to the optimal barcode set with distance", optimal.dist))
        } else {
		optimal.barcodes <- kept.barcodes
		optimal.dist <- mean.dist
		print(paste("the current average distance of", mean.dist, "is higher than the optimal barcode set with distance", optimal.dist))
        }
}

if (iter == iter.number) {
	# if number of desired barcodes was not met
	if (length(optimal.barcodes) == 0) {
		stop(paste("desired barcode number of", desired.barcodes, "could not be achieved"))
	} else {
	# if number of desired barcodes was met
		print(".....exporting optimal barcode set.....")
		print(optimal.barcodes)
	
		barcodes.export <- tibble(
			"Barcode_ID" = sprintf("BC%04d", 1:desired.barcodes),
			"Barcodes" = optimal.barcodes,
			"Barcode_Reverse_Complement" = map_chr(optimal.barcodes, ~reverse.complement(.)),
			"Adapter_Sequence" = map_chr(
				optimal.barcodes, 
				~ paste0(
					"CAAGCAGAAGACGGCATACGAGAT",
					reverse.complement(.),
					"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNN"
				)
			)
		)
	
		write.table(
			barcodes.export, "barcodes.tsv",
			sep="\t", col.names=T, row.names=F, quote=F, na=""
		)
	
		print(paste(optimal.dist, " average hamming distance for optimal barcode set"))
		stop(paste("desired barcode number", desired.barcodes, "met, output saved"))
	}
}

}
