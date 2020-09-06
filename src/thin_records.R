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


# Load Data ---------------------------------------------------------------
wtjr.occ.unique.gps <- read.csv(args[1], stringsAsFactors = F)

# Do thinning -------------------------------------------------------------
set.seed(496)
wtjr.occ.unique.gps %>% 
  mutate(SPEC = "Leptown") %>% 
  thin(., thin.par = 1,
       lat.col = "roundlat", long.col = "roundlon", 
       reps = 3*dim(.)[1], out.dir = "processed_data/thin/1km/",
       out.base = "wtjr_occ", log.file = "processed_data/thin/1km/thin_log.txt")
wtjr.occ.unique.gps %>% 
  mutate(SPEC = "Leptown") %>% 
  thin(., thin.par = 5,
       lat.col = "roundlat", long.col = "roundlon", 
       reps = 3*dim(.)[1], out.dir = "processed_data/thin/5km/",
       out.base = "wtjr_occ", log.file = "processed_data/thin/5km/thin_log.txt")
wtjr.occ.unique.gps %>% 
  mutate(SPEC = "Leptown") %>% 
  thin(., thin.par = 10,
       lat.col = "roundlat", long.col = "roundlon", 
       reps = 3*dim(.)[1], out.dir = "processed_data/thin/10km/",
       out.base = "wtjr_occ", log.file = "processed_data/thin/10km/thin_log.txt")
wtjr.occ.unique.gps %>% 
  mutate(SPEC = "Leptown") %>% 
  thin(., thin.par = 50,
       lat.col = "roundlat", long.col = "roundlon", 
       reps = 3*dim(.)[1], out.dir = "processed_data/thin/50km/",
       out.base = "wtjr_occ", log.file = "processed_data/thin/50km/thin_log.txt")


