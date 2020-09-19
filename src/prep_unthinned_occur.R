####################
# Project: WTJR
# Author: Timothy Thurman
# Purpose: Make thinned datasets for SDMs
# Date Created: Fri Sept  4 20:05:59 2020
####################

# Meant to be called from the Snakemake pipeline

# Load packages -----------------------------------------------------------
library(tidyverse)


# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 input file
# 2 output file


# Load Data ---------------------------------------------------------------
occur <- read.csv(args[1], stringsAsFactors = F) %>%
    mutate(SPEC = "Leptown") %>%
    select(SPEC, roundlon, roundlat)
write.csv(occur, file = args[2], row.names = F, quote = F)