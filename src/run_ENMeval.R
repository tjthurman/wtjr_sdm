####################
# Project: WTJR
# Author: Timothy Thurman
# Purpose: Running SDMs with ENMeval
# Date Created: Fri Apr 24 12:39:58 2020
####################


# Load packages -----------------------------------------------------------
library(tidyverse)
library(raster)
library(dismo)
library(ENMeval)


# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 csv of occurrence records for analysis
# 2 bioclim data
# 3 bg points
# 3 feature class
# 4 seed
# 5 cores
in_csv <- args[1]
in_bioclim <- args[2]
bg <- args[3]
fc <- args[4]
seed <- args[5]
cores <- args[6]

# Load occurence data -----------------------------------------------------
wtjr.occ <- read.csv(in_csv, stringsAsFactors = F) %>%
  dplyr::select(roundlon, roundlat)

# Load Bioclim data -------------------------------------------------------
bioclim <- brick(in_bioclim)
# Fix names
names(bioclim) <- paste("bio", 1:19, sep = "")
extent <- extent(bioclim)

# Background Points -----------------------------------------------------------------
load(bg)

# Parse input filename ----------------------------------------------------
dataset <- str_remove(in_csv, pattern = ".csv") %>% 
  str_split(., pattern = "_") %>% 
  .[[1]] %>%
  .[4:5] %>% 
  paste(., sep = "", collapse = "_")
  
# Loop to allow multiple running attempts ---------------------------------
set.seed(seed)

# Some datasets get consistent convergence errors in GLMNET,
# which runs the model fitting in ENMeval. 

# I'll try using a different thinned dataset:
# I used the first thinned dataset by default. Snakemake is bad for this,
# so instead I'll just do the switch in here.

# For the two thinning distances that are having convergence troubles in some 
# feature classes (1km and 5km), I'll use the second
# .csv output in the thinning, instead of the first. 
if (dataset %in% c("5km_thin1")) {
  in_csv2 <- str_replace(in_csv, pattern = "thin1", "thin2")
  
  set.seed(245)
  wtjr.occ <- read.csv(in_csv2, stringsAsFactors = F) %>%
    dplyr::select(roundlon, roundlat) 
  print("Using csv file:")
  print(in_csv2)
}

if (dataset %in% c("1km_thin1")) {
  in_csv2 <- str_replace(in_csv, pattern = "thin1", "thin3")
  
  set.seed(245)
  wtjr.occ <- read.csv(in_csv2, stringsAsFactors = F) %>%
    dplyr::select(roundlon, roundlat) 
  print("Using csv file:")
  print(in_csv2)
}

attempt_limit <- 5
target.file <- paste("results/enmeval_res_", dataset, "_", fc, ".RData", sep = "")

attempt <- 1
while (attempt <= attempt_limit) {
  # move the attempt counter up
  attempt <- attempt + 1
  if (!file.exists(target.file)) { # If the results file doesn't exist
    # Run the model within a try statement and save the results
    print(paste0("Starting ENM eval attempt ", attempt - 1))
    enmeval_res <- try(ENMevaluate(occ = wtjr.occ,
                                   env = bioclim,
                                   bg.coords = background,
                                   algorithm = "maxnet",
                                   method = "checkerboard2",
                                   parallel = F, rasterPreds = T, 
                                   numCores = NULL, updateProgress = T,
                                   fc = fc))
    print(paste0("Finished ENM eval attempt ", attempt - 1))
    print("start saving")
    save(enmeval_res, file = target.file)
    print("end saving")
    
    # Then check to see if the result file is big enough
    if (file.size(target.file) < 1e6) { # if not
      # remove it and increment the attempt counter
      print(paste0("ENM eval attempt ", attempt - 1, " failed, retrying"))
      file.remove(target.file)
    } else {
      print(paste0("ENM eval attempt ", attempt - 1, " worked"))
      break
    }
  } else {
    break
  }
}
