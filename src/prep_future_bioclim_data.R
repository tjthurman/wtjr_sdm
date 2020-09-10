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


# Do  processing ------------------------------------------------------------
# Loop over bioclim variables (need 2 and 3)
# Then loop over folders within the base directory

for (bioclim_var in paste0("bio_", 2:3)) { # For each bioclim variable
  x <- list() # set up a list of results
  i <- 1
  for (folder in list.files(base_dir)) { # loop over the 5 folders containing the projections from the 5 models
    layer <- raster(x = file.path(base_dir, folder, bioclim_var)) # grab the raster of the right bioclim layer
    x[[i]] <- layer # add it to the list
    i <- i + 1
  }
  all_models_cropped <- crop(stack(x), y = template) # then convert list to a raster stack and crop to right extent
  model_avg <- calc(x = all_models_cropped, fun = mean, na.rm = T) # take the model averaged mean
  assign(x = eval(bioclim_var), value = model_avg) # and assign it to an object named according to the bioclim variable
} 

# COmbine into rasterstack
stacked <- stack(bio_2, bio_3)

# And save
writeRaster(x = stacked, filename = output_file)

