# Code and Data for

## XXX TITLE OF MS TO COME

### About

This Github repository contains the code, and some of the data, needed to recreate the species distribution modelling (SDM) and phenotypic modelling analyses of:


XXX CITATION FOR PAPER HERE. 

We have used a set of software tools to try to make this analysis pipeline as reproducible as possible. We use the Python program `Snakemake` to create a reproducible analysis pipeline. Most analyses are performed in `R`, and we use the `renv` package to manage libraries/packages for `R`. 

### Installation

You will need the following software installed on your system:

1. Python 3. We recommend installing Python 3 through [Anaconda](https://www.anaconda.com) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html). We use Python version XXX LIST IT HERE. 
2. Snakemake. We recommend following instaling `snakemake` through `conda`, see [these instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 
3. R. R can be download from the [R Project for Statisical Computing](https://www.r-project.org). We used R version `4.0.2`. 
4. The R package `renv`. See [here](https://rstudio.github.io/renv/index.html) for instructions on installing `renv`. 5. R package dependencies. These dependencies are managed with the `renv` package. To install the dependencies, open this project in RStudio (easiest) or, in your terminal, set this project as your working directory and then launch `R`. Then, call `renv::restore()` to install the necessary packages (this may take a while). 

### Run the pipeline

Once everything is installed, you can use snakemake to run the analysis pipeline by simply executing:

```
snakemake --cores 1 
```

Note that this will take a long time. 

### Structure

This repository is set up as an RStudio project. You can open it in RStudio by opening the `wtjr_sdm.Rproj` file. The project is organized into a a few top-level folders:


* `docs`- Many of the analysis scripts in our pipeline are `.Rmd` documents that mix code and text. These are rendered into `.html` documents with  `rmarkdown` and `knitr`. The `docs` folder contains the output `.html` document corresponding to the `.Rmd` file of the same name in the `src` directory. 
* `processed_data`- Data that has been processed/curated and is ready for analysis. All these files are created from raw data file (in the `raw_data` folder) and analysis scripts found in `src/`.
* `raw_data`- Raw data obtained for this project, without any processing. See README in that folder for more details on data source for each item. 
* `renv`- package/library management for this project. This folder is created/maintained by the `renv` package, users should not need to edit it. 
* `results`- Analysis results, generated from raw data, processed data, and scripts. 
8. `src`- Scripts and functions written for this project. Some script are in `.Rmd` files, which generate corresponding `.html` files in the `docs` directory. 



Old files


2. bin- Code (and some data) taken from other WTJR projects. The `Ecological_modelling_outputs/` folder contains Alex Kumar's work on phenotypic modelling re: conservaion status. The `mills_2018_scripts_and_data` folder contains the data from the Mills et al/ 2018 paper on winter white coloration, downloaded from this Dryad link: https://datadryad.org/stash/dataset/doi:10.5061/dryad.8m0p1. However, I have deleted and not tracked the `rasterstack4cov.tif` file of global environmental predictors, as it is very large. A version of this file, cropped down to just what we need for WTJR, is available in `processed_data/`.


