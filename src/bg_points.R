# Load packages -----------------------------------------------------------
library(dismo)
library(raster)

# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 bioclim data
# 2 seed
in_bioclim <- args[1]
seed <- args[2]

# Load Bioclim data -------------------------------------------------------
bioclim <- brick(in_bioclim)
# Fix names
names(bioclim) <- paste("bio", 1:19, sep = "")
extent <- extent(bioclim)

# Background Points -----------------------------------------------------------------
# Generate a common set of background points
set.seed(seed)
background <- randomPoints(mask = bioclim,     
                           n = 10000, 
                           extf = 1)
save(background, file = "processed_data/bg_points_for_sdm.RData")
