####################
# Project: WTJR
# Author: Timothy Thurman
# Purpose: Predict current color phenotype from environmental data
# Date Created: Thu Sept 25 2:46:58 2020
####################


# Load packages -----------------------------------------------------------
library(tidyverse)
library(readxl)
library(raster)



# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 phenotype data
# 2 environment data
# 3 gbif data
# 4 SDM range

pheno_dat <- args[1]
environ <- args[2]
gbif_dat <- args[3]
range_file <- args[4]

# Load data -----------------------------------------------------------
# Pheno data
pheno <- read_excel(pheno_dat)  %>% 
  mutate(colorSymbol = ifelse(`binary color` == 2, 0, 1)) %>% 
  rename(snow.old = snow, bio2_old = bio_2, bio3_old = bio_3, alt_old = alt)

# Environmental predictors
predictors <- brick(environ)
names(predictors)<-c("snow.cover", "bio_2", "bio_3")

# Load info from GBIF
gbif <- read.delim(gbif_dat, quote = "", stringsAsFactors = F) %>% 
  filter(catalogNumber %in% c(2869, 2870))

rangemap <- raster(range_file)

# Analysis -----------------------------------------------------------

# Fix location info for the two incorrect samples
pheno$decLat[str_detect(pheno$source, "MCZ")] <- as.character(gbif$decimalLatitude)
pheno$decLong[str_detect(pheno$source, "MCZ")] <- as.character(gbif$decimalLongitude)

# Extract environmental data for locations with phenotuped individuals
env <- raster::extract(predictors, y = SpatialPoints(data.frame(x = as.numeric(pheno$decLong),
                                                                y = as.numeric(pheno$decLat)), 
                                                                proj4string = CRS(proj4string(predictors))))
# Add it to our pheno data
pheno <- cbind(pheno, env) 

# Run
lepto.mod <- glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = "binomial", data = pheno, na.action = na.pass) 
# Make predictions
lepto.mod.pred <- predict(predictors, lepto.mod, type='response')

# Mask predictions to range
pheno.range <- mask(lepto.mod.pred, rangemap)

# Save results ------------------------------------------------------------
save(lepto.mod, file = "results/pheno/current_pheno_glm.RData")
writeRaster(lepto.mod.pred, filename = "results/pheno/current_predicted_probWhite", format='GTiff')

writeRaster(pheno.range, filename = "results/pheno/current_predicted_probWhite_SDMrange", format='GTiff')
