# Barcode Generator

## About

This software will let you generate barcoded Reverse Transcription Oligos (RTOs) for STRIPE-seq that are optimized for Illumina short read sequencing.

## Getting Started

### Cloning Repository

To get started, you must first clone the Barcode Generator repository.
Navigate to a directory you would like to clone the repository to and enter `https://github.com/rpolicastro/barcode-generator.git`.

### Prepare Software

The Barcode Generator package takes advantage of the R libraries tidyverse and gtools.
If you already have R installed and these packages present, you are ready to run the software.
If not, it is recommended to install the software either through the conda software package manager, or to use the provided singularity container.

#### Conda Environment

This software was originally developed using the [conda](https://conda.io/en/latest/) package manager and virtual environment.
The conda package manager installs both the main software and all dependencies into a 'virtual environment' to ensure compatabilty.

Before creating the environment, you must first install miniconda.
1. [Install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda), and make sure that conda is in your PATH.
2. Update conda to the latest version `conda update -n base -c defaults conda`.

After miniconda is installed, you are now ready to create the conda environment with the appropriate software.

1. Create the new environment and specify the software to include in it.
```
conda create -n barcode-generator -y -c conda-forge -c bioconda \
r-tidyverse r-gtools
```
2. Update the software to the latest compatible versions.
```
conda update -n barcode-generator -y -c conda-forge -c bioconda --all
```

To use the software in the environment you can type `conda activate barcode-generator`. You can deactivate the environment by closing your terminal or entering `conda deactivate`.

#### Singularity Container

Singularity containers are self contained 'boxes' that houses the software required to run the Barcode Generator script. Before running the script you must install the Singularity software, and download the Barcode Generator container.

1. Install the [latest version](https://sylabs.io/guides/3.3/user-guide/installation.html) of Singularity.
2. Grab the Barcode Generator singularity container.
```
singularity pull library://rpolicastro/default/barcode_generator:1.0.0
```

To use the software within the container, you must enter the container and specify which directories you will be working in, as well as load the interncal conda environment.
1. Shell into the container.
```
singularity shell \
-e -C \
-B path/to/repository \
-H path/to/repository \
path/to/container/barcode_generator_1.0.0.sif
```
2. Activate the internal conda environment `source activate barcode-generator`

# Built With

This workflow would not be possible without the great software listed below.

- [Anaconda](https://www.anaconda.com/) - Software package manager.
- [R](https://www.r-project.org/) - Robust language for statistical computing.
- [Tidyverse](https://www.tidyverse.org/) - Data manipulation in R.
- [gtools](https://cran.r-project.org/web/packages/gtools/index.html) - R library with various conveneince functions.
