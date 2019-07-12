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
conda update -n chipseq-automation -y -c conda-forge -c bioconda --all
```

To use the software in the environment you can type `conda activate barcode-generator`. You can deactivate the environment by closing your terminal or entering `conda deactivate`.
