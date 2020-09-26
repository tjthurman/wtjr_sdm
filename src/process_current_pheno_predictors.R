####################
# Project: WTJR
# Author: Timothy Thurman
# Purpose: Crop the predictors from Mills et al. 2018 into something smaller and more useful for our project. 
# Date Created: Tue Apr 21 15:05:37 2020
####################


# Load packages -----------------------------------------------------------
library(tidyverse)
library(raster)

# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 input snow file
# 2 input bioclim file
# 3 output

mills <- args[1]
bc <- args[2]
output <- args[3]

# Read in predictors file from mills -------------------------------------------------
mills_preds <- stack(mills)
names(mills_preds) <- c("snow.cover", "bio_2", "bio_3","alt")
snow <- subset(mills_preds, subset = "snow.cover")
crop.snow <- raster::crop(x = snow, y = extent(c(-130, -85, 30, 65)))



# Read in bioclim -------------------------------------------------
bioclim <- stack(bc)
names(bioclim) <- paste0("bio_", 1:19)
bc23 <- subset(bioclim, c("bio_2", "bio_3"))

# Combine
pheno.preds <- stack(crop.snow, bc23)


# Crop -----------------------------------------------
crop.pheno.preds <- raster::crop(x = pheno.preds, y = extent(c(-130, -85, 30, 65)))

# Save as multi-band tiff -------------------------------------------------
writeRaster(crop.pheno.preds, filename=output,
            format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
