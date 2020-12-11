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
# 1 input phenotype data
# 2 input environment data
# 3 input gbif data
# 4 input SDM range
# 5 output RData of glm
# 6 output probwhite across full raster
# 7 output probwhite across SDM range

pheno_dat <- args[1]
environ <- args[2]
gbif_dat <- args[3]
range_file <- args[4]
out_glm <- args[5]
out_raster_full <- args[6]
out_raster_range <- args[7]

# For running as script
# pheno_dat <- "raw_data/Ltowsendii_database_FINAL.xlsx"
# environ <- "processed_data/pheno_predictors_millsetal2018.tif"
# gbif_dat <- "raw_data/GBIF/verbatim.txt"
# range_file <- "results/sdm/sdm_rangemap_best_sens95.grd"
# out_glm <- "results/pheno/current_pheno_glm.RData"
# out_raster_full <- "results/pheno/current_predicted_probWhite.tif"
# out_raster_range <- "results/pheno/current_predicted_probWhite_SDMrange.tif"


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

plot(lepto.mod.pred)
plot(pheno.range)
# Save results ------------------------------------------------------------
save(lepto.mod, file = out_glm)
writeRaster(lepto.mod.pred, filename = out_raster_full, format='GTiff', overwrite = T)

writeRaster(pheno.range, filename = out_raster_range, format='GTiff', overwrite = T)
