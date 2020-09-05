####################
# Project: WTJR
# Author: Timothy Thurman
# Purpose: prepare Bioclim data for WTJR SDM
# Date Created: Tue Apr  7 14:31:59 2020
####################

# Load packages -----------------------------------------------------------
library(tidyverse)
library(raster)
library(rgdal)
library(stringr)



# Download worldclim data -------------------------------------------------
# We'd like to have data for the extent from:
# West= -130, 
# East = -85, 
# South = 30
# North = 65

# Unfortunately, this extent is spread across 6 tiles. 
# So, we want to download all the tiles, crop them down to the area we want, and then merge them.

# North row, west to east
tile.01 <- getData(name = "worldclim",
                   var = "bio",
                   res = 0.5,
                   lon = -130,
                   lat = 64,
                   path = "raw_data/")
tile.02 <- getData(name = "worldclim",
                   var = "bio",
                   res = 0.5,
                   lon = -110,
                   lat = 64,
                   path = "raw_data/")
tile.03 <- getData(name = "worldclim",
                   var = "bio",
                   res = 0.5,
                   lon = -85,
                   lat = 64,
                   path = "raw_data/")

extent(tile.01)
extent(tile.02)
extent(tile.03)


# Southern row, west to east
tile.11 <- getData(name = "worldclim",
                           var = "bio",
                           res = 0.5,
                           lon = -130,
                           lat = 46,
                           path = "raw_data/")
tile.12 <- getData(name = "worldclim",
                             var = "bio",
                             res = 0.5,
                             lon = -112,
                             lat = 46,
                             path = "raw_data/")
tile.13 <- getData(name = "worldclim",
                   var = "bio",
                   res = 0.5,
                   lon = -85,
                   lat = 46,
                   path = "raw_data/")


extent(tile.11)
extent(tile.12)
extent(tile.13)


# Crop tiles to desired extent --------------------------------------------
tile.01.crop <- crop(tile.01, y = extent(x = c(-130, -120, 60, 65)))
tile.02.crop <- crop(tile.02, y = extent(x = c(-120, -90, 60, 65)))
tile.03.crop <- crop(tile.03, y = extent(x = c(-90, -85, 60, 65)))

tile.11.crop <- crop(tile.11, y = extent(x = c(-130, -120, 30, 60)))
tile.12.crop <- tile.12
tile.13.crop <- crop(tile.13, y = extent(x = c(-90, -85, 30, 60)))

# Merge cropped tiles -----------------------------------------------------
cropped.tiles <- list(tile.01.crop,
                      tile.02.crop,
                      tile.03.crop,
                      tile.11.crop,
                      tile.12.crop,
                      tile.13.crop)
final.layer <- do.call(merge, cropped.tiles)


# Fix names ---------------------------------------------------------------
names(final.layer)
layer.names <- str_remove(names(tile.01.crop), pattern = "\\_01")
names(final.layer) <- layer.names


# Save as multi-band tiff -------------------------------------------------
writeRaster(final.layer, filename="processed_data/bioclim_30arcsec_for_WTJR_SDM.tif",
            format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
