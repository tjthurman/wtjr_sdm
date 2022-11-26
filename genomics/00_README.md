This folder contains code and data needed to recreate genomic analyses from:

XXX CITATION FOR PAPER HERE. 

## Folder organization

This folder is organized according to the analysis steps we performed:

* `01_process_raw_reads/` - Code for processing raw read data this code must be run before all other analyses. See README within for more information on running the analysis and software environments. 

* `02_genome-wide_PCA/`- Code for performing a genome-wide PCA and creating the PCA plot presented in the supplemental material. See README within for more information on running the analysis and software environments. 

* `03_color-polymoprhism-across-range/`- Code for estimating allele frequencies of the top color-associated SNPs inside and outside Colorado and generating the piechart presented in the supplemental material. See README within for more information on running the analysis and software environments. 

* `04_Fst_CO-vs-NDK/` - Code for estimating genome-wide Fst between white-tailed jackrabbit samples from Colorado and North Dakota. See README within for more information on running the analysis and software environments. 

* `05_dosage-environment-correlation/`- Code for estimating the correlation between dosage of white-associated allele and environmental variables across the range. See README within for more information on running the analysis and required software. 


## Data dependencies

All raw sequencing data for these analyses are available from the NCBI short-read archive at under BioProject number PRJNA726805: https://www.ncbi.nlm.nih.gov/sra/PRJNA726805. To replicate our analyses, users may need to change filepaths in our code to be compatible with their file system. 


## Software 

We used snakemake to manage workflows for these analyses, and used `conda` and `renv` to manage software dependencies. For each analysis, we include a `conda` environment file describing the software dependencies needed. These analyses also require R (we used version 4.0.2). See top-level repository for more information on installing these software packages. 


