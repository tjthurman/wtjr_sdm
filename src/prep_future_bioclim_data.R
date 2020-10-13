####################
# Project: WTJR
# Author: Timothy Thurman
# Purpose: Prepping projected Bioclim data
# to match the SRT data from the forest service
# Date Created: Wed Sep 9 3:59:36 2020
####################


# Load packages -----------------------------------------------------------
library(raster)


# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 directory containing raw CMIP5 data
base_dir <- args[1]
# 2 current bioclim data to be used as extent object for cropping
template <- brick(args[2])
# 3 output file
output_file <- args[3]

base_dir <- "raw_data/worldclim_projections_2080s/"
template <- brick("processed_data/bioclim_30arcsec_for_WTJR_SDM.tif")
output_file <- "processed_data/bc23_CMIP5_RCP85_2080s_5modavg.grd"

# Do  processing ------------------------------------------------------------
# Loop over bioclim variables (need 2 and 3)
# Then loop over folders within the base directory
for (bioclim_var in paste0("bio_", 2:3)) { # For each bioclim variable
  x <- list() # set up a list of results
  i <- 1
  for (folder in list.files(base_dir)) { # loop over the 5 folders containing the projections from the 5 models
    layer <- raster(x = file.path(base_dir, folder, paste0(bioclim_var, ".asc"))) # grab the raster of the right bioclim layer
    crop_layer <- crop(layer, y = template)
    x[[i]] <- crop_layer # add it to the list
    i <- i + 1
    rm(layer, crop_layer)
  }
  all_models_cropped <- stack(x) # then convert list to a raster stack
  model_avg <- calc(x = all_models_cropped, fun = mean, na.rm = T) # take the model averaged mean
  assign(x = eval(bioclim_var), value = model_avg) # and assign it to an object named according to the bioclim variable
} 

# Combine into rasterstack
stacked <- stack(bio_2, bio_3)

# And save
writeRaster(x = stacked, filename = output_file, overwrite = T)

