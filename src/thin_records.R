####################
# Project: WTJR
# Author: Timothy Thurman
# Purpose: Make thinned datasets for SDMs
# Date Created: Fri Sept  4 20:05:59 2020
####################

# Meant to be called from the Snakemake pipeline

# Load packages -----------------------------------------------------------
library(tidyverse)
library(spThin)



# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 input file
# 2 thinning distance
# 3 reps
# 4 out directory
# 5 seed


# Load Data ---------------------------------------------------------------
wtjr.occ.unique.gps <- read.csv(args[1], stringsAsFactors = F)

# Do thinning -------------------------------------------------------------
set.seed(as.integer(args[5]))
wtjr.occ.unique.gps %>% 
  mutate(SPEC = "Leptown") %>% 
  thin(., thin.par = as.numeric(args[2]),
       lat.col = "roundlat", long.col = "roundlon", 
       reps = as.numeric(args[3]), out.dir = args[4],
       out.base = paste0("wtjr_occ_", args[2], "km"), log.file = file.path(args[4], paste0("thin_log_", args[2], "km.txt")))









