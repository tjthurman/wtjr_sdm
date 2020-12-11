####################
# Project: WTJR
# Author: Timothy Thurman
# Purpose: Predict future color phenotype from environmental data
# Date Created: Thu Sept 25 2:46:58 2020
####################


# Load packages -----------------------------------------------------------
library(tidyverse)
library(readxl)
library(raster)



# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 phenotype data
# 2 current bc
# 3 current srt
# 4 future bc
# 5 future srt
# 6 SDM range
# 7 gbif
# 8 output glm
# 9 output current SRT raster
# 10 output future SRT raster
pheno_dat <- args[1]
cur_bc_file <- args[2]
cur_srt_file <- args[3]
future_bc_file <- args[4]
future_srt_file <- args[5]
range_file <- args[6]
gbif_dat <- args[7]
glm_out <- args[8]
current_SRT <- args[9]
future_SRT <- args[10]

# For running as script
# pheno_dat <- "raw_data/Ltowsendii_database_FINAL.xlsx"
# cur_bc_file <- "processed_data/bioclim_30arcsec_for_WTJR_SDM.tif"
# cur_srt_file <- "raw_data/SRT/SRT_historical/SRT_historical.tif"
# future_bc_file <- "processed_data/bc23_CMIP5_RCP85_2080s_5modavg.grd"
# future_srt_file <- "raw_data/SRT/SRT_RCP85_2080s/SRT_RCP85_2080s.tif"
# range_file <- "results/sdm/sdm_rangemap_best_sens95.grd"
# gbif_dat <- "raw_data/GBIF/verbatim.txt"
# glm_out <- "results/pheno/current_pheno_glm_SRT.RData"
# current_SRT <- "results/pheno/current_predicted_probWhite_SDMrange_SRT.tif"
# future_SRT <- "results/pheno/future_predicted_probWhite_SDMrange.tif"

# Load data -----------------------------------------------------------
# Pheno data
pheno <- read_excel(pheno_dat)  %>% 
  mutate(colorSymbol = ifelse(`binary color` == 2, 0, 1)) %>% 
  rename(snow.old = snow, bio2_old = bio_2, bio3_old = bio_3, alt_old = alt)

# Load info from GBIF
gbif <- read.delim(gbif_dat, quote = "", stringsAsFactors = F) %>% 
  filter(catalogNumber %in% c(2869, 2870))

# Fix location info for the two incorrect samples
pheno$decLat[str_detect(pheno$source, "MCZ")] <- as.character(gbif$decimalLatitude)
pheno$decLong[str_detect(pheno$source, "MCZ")] <- as.character(gbif$decimalLongitude)


# Environmental predictors
current.srt <- raster(cur_srt_file)
future.srt <- raster(future_srt_file)
current.bc <- brick(cur_bc_file)
names(current.bc) <- paste0("bio_", 1:19)
current.bc23 <- raster::subset(current.bc, subset = c("bio_2", "bio_3"))
future.bc23 <- brick(future_bc_file)
names(future.bc23) <- paste0("bio_", 2:3)

# Disaggregate SRT to high resolution
current.srt.disagg <- disaggregate(current.srt, fact = 5)
future.srt.disagg <- disaggregate(future.srt, fact = 5)

# Range map
range <- raster(range_file)




# Re-do phenotypic model using current SRT ----------------------------------------
# Extract the SRT info for each location
srt <- raster::extract(current.srt.disagg, y = SpatialPoints(data.frame(x = as.numeric(pheno$decLong),
                                                                 y = as.numeric(pheno$decLat)), 
                                                      proj4string = CRS(proj4string(current.srt.disagg))))
# Extract the bioclim variables for each location
bc23 <- raster::extract(current.bc23, y = SpatialPoints(data.frame(x = as.numeric(pheno$decLong),
                                                                 y = as.numeric(pheno$decLat)), 
                                                      proj4string = CRS(proj4string(current.bc23))))
# Combine together for the GLM
glm.data.current <- cbind(pheno, srt, bc23) %>% 
  filter(!is.na(srt))
srt.mod <- glm(colorSymbol ~ srt + bio_2 + bio_3, family = "binomial", data = glm.data.current, na.action = na.pass) 

# Predict current probWhite with SRT --------------------------------------
# Combine current rasters
tmp1 <- raster::intersect(current.srt.disagg, current.bc23)
tmp2 <- raster::intersect(current.bc23, current.srt.disagg)
extent(tmp1) <- alignExtent(tmp1, tmp2)
current.env <- brick(tmp1, tmp2)
names(current.env) <- c("srt", "bio_2", "bio_3")

current.pred.pheno.srt <- predict(current.env, srt.mod, type = "response")
current.pred.pheno.srt.range <- mask(current.pred.pheno.srt, crop(range, current.pred.pheno.srt))


# Predict future probWhite with SRT ------------------------------------------------
tmp3 <- raster::intersect(future.srt.disagg, future.bc23)
tmp4 <- raster::intersect(future.bc23, future.srt.disagg)
extent(tmp3) <- alignExtent(tmp3, tmp4)
future.env <- brick(tmp3, tmp4)
names(future.env) <- c("srt", "bio_2", "bio_3")


future.pred.pheno <- predict(future.env, srt.mod, type = "response")
future.pred.pheno.range <- mask(future.pred.pheno, crop(range, future.pred.pheno))


# Save results ------------------------------------------------------------
# Save the SRT model
save(srt.mod, file = glm_out)

# Save the current phenotypic predictions using SRT
writeRaster(current.pred.pheno.srt.range, filename = current_SRT, format = "GTiff", overwrite = T)

# Save the future predictions using SRT
writeRaster(future.pred.pheno.range, filename = future_SRT, format = "GTiff", overwrite = T)
