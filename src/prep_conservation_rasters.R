####################
# Project:
# Author: Timothy Thurman
# Purpose:
# Date Created: Tue Apr 21 15:05:30 2020
####################

# Load packages -----------------------------------------------------------
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(sp)
library(raster)

# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 phenotype data
# 2 environment data
# 3 gbif data
# 4 SDM range

bioclim <- args[1]



# Get Shapefiles of states and provinces ----------------------------------
state_prov <- ne_states(c("united states of america", "canada"))

# Set up name lists -------------------------------------------------------
# according to Fig. 5 of Brown at al. 2020
extirpated.states <- c("Oregon", "British Columbia", "Kansas","Oklahoma", "Missouri", "Illinois", "Wisconsin")
broad.extirp.states <- c("California", "Washington", "Nebraska", "Iowa", "Minnesota")
local.extirp.states <- c("Nevada", "Utah", "Colorado", "New Mexico", "Wyoming", "South Dakota")
poss.decline.states <- c("Idaho", "Montana", "Alberta", "Saskatchewan", "North Dakota")
pres.stable.stats <- c("Manitoba")


# Subset down to polygons -------------------------------------------------
extirpated.poly <- state_prov[state_prov$name %in% extirpated.states,]
broad.extirp.poly <- state_prov[state_prov$name %in% broad.extirp.states,]
local.extirp.poly <- state_prov[state_prov$name %in% local.extirp.states,]
poss.decline.poly <- state_prov[state_prov$name %in% poss.decline.states,]
pres.stable.poly <- state_prov[state_prov$name %in% pres.stable.stats,]


# Set up a template -------------------------------------------------------
# Use the bioclim data as a template, since it has the full extent of what we're looking at
temp <- brick(bioclim)



# Make the rasters --------------------------------------------------------
extirpated.raster <- rasterize(extirpated.poly, temp, field = 1)
broad.extirp.raster <- rasterize(broad.extirp.poly, temp, field = 1)
local.extirp.raster <- rasterize(local.extirp.poly, temp, field = 1)
poss.decline.raster <- rasterize(poss.decline.poly, temp, field = 1)
pres.stable.raster <- rasterize(pres.stable.poly, temp, field = 1)


# Save the rasters --------------------------------------------------------
writeRaster(extirpated.raster, filename="processed_data/conservation_rasters/extirpated.tif",
            format="GTiff", overwrite=TRUE,options=c("COMPRESS=LZW"))
writeRaster(broad.extirp.raster, filename="processed_data/conservation_rasters/broad_extirp.tif",
            format="GTiff", overwrite=TRUE,options=c("COMPRESS=LZW"))
writeRaster(local.extirp.raster, filename="processed_data/conservation_rasters/local_extirp.tif",
            format="GTiff", overwrite=TRUE,options=c("COMPRESS=LZW"))
writeRaster(poss.decline.raster, filename="processed_data/conservation_rasters/poss_decline.tif",
            format="GTiff", overwrite=TRUE,options=c("COMPRESS=LZW"))
writeRaster(pres.stable.raster, filename="processed_data/conservation_rasters/pres_stable.tif",
            format="GTiff", overwrite=TRUE,options=c("COMPRESS=LZW"))


