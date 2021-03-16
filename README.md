Code and Data for

## XXX TITLE OF MS TO COME

### About

This Github repository contains the code, and some of the data, needed to recreate the species distribution modelling (SDM) and phenotypic modelling analyses of:

XXX CITATION FOR PAPER HERE. 

We have used a set of software tools to try to make this analysis pipeline as reproducible as possible. We use the Python program `Snakemake` to create a reproducible analysis pipeline. Most analyses are performed in `R`, and we use the `renv` package to manage libraries/packages for `R`. 

### Installation

You will need the following software installed on your system:

1. Python 3. We recommend installing Python 3 through [Anaconda](https://www.anaconda.com) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html). We used Python version 3.8.3. 
2. Snakemake. We recommend following installing `snakemake` through `conda`, see [these instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). We used version 5.19.3.
3. R. R can be downloaded from the [R Project for Statisical Computing](https://www.r-project.org). We used R version `4.0.2`. 
4. The R package `renv`. See [here](https://rstudio.github.io/renv/index.html) for instructions on installing `renv`. 5. R package dependencies. These dependencies are managed with the `renv` package. To install the dependencies, open this project in RStudio (easiest) or, in your terminal, set this project as your working directory and then launch `R`. Then, call `renv::restore()` to install the necessary packages (this may take a while). 

### Structure

This repository is set up as an RStudio project. You can open it in RStudio by opening the `wtjr_sdm.Rproj` file. The project is organized into a  few top-level folders:

* `processed_data`- Data that has been processed/curated and is ready for analysis. All these files are created from raw data file (in the `raw_data` folder) and analysis scripts found in `src/`.
* `raw_data`- Raw data obtained for this project, without any processing. See the README file in that folder for more details on the data source for each item. **NB**- Many raw data files are included in this repository, but some raw data files are too large to be uploaded to Github. See the README file in the raw data folder for information on how to obtain those files, either from Dryad or from the original sources we used. 
* `renv`- package/library management for this project. This folder is created/maintained by the `renv` package, users should not need to edit it. 
* `results`- Analysis results, generated from raw data, processed data, and scripts. 
* `src`- Scripts and functions written for this project. Some script are in `.Rmd` files, which generate corresponding `.html` files upon knitting. 

### Data dependencies

Some raw data files are too large to be uploaded to Github. They muse be obtained through other sources and then copied into the appropriate folder before running the pipeline. See the README file in the `raw_data/` folder for instructions on how to obtain these files. 

### Run the pipeline

Once all the necessary software is installed and you have obtained the necessary data, you can use `snakemake` to run our analysis pipeline. To do this, set this project directoy as your working directory, and then execute:

```
snakemake --cores 1 
```

This should run the full pipeline, though note that it will likely take a long time, and may fail depending on the amount of procesing power and memory available on your computer. We ran our analysis on a HPC cluster, and recommend you do the same. See the `snakemake` documentation for guidance on running `snakemake` on a [cluster](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) and how to set up a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) for specifying HPC options. A sample command for running our pipeline on a cluster, using the `snakemake` profile found in a folder named `sp_wtjr`, is:

```
snakemake --profile sp_wtjr/ 
```

but you will need to create your own profile adapted to your own HPC. 

