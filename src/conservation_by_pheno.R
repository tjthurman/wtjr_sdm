####################
# Project: WTJR SDMs
# Author: Timothy Thurman
# Purpose: analysis of phenotype by conservation status overlap
# Date Created: Tue Apr 21 15:05:30 2020
####################

# Load packages -----------------------------------------------------------
library(tidyverse)
library(raster)
library(broom)

# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 input phenotype data
# 2 output file
pheno.file <- args[1]
outfile <- args[2]

# For running as script
# pheno.file <- "results/pheno/current_predicted_probWhite_SDMrange.tif"
# outfile <- "results/conservation/cons_by_current_color.RData"


# Load pheno raster -----------------------------------------------
pheno.range <- raster(pheno.file)

# Load conservation rasters -----------------------------------------------
extirpated <- raster("processed_data/conservation_rasters/extirpated.tif")
broad.extirp <- raster("processed_data/conservation_rasters/broad_extirp.tif")
local.extirp <- raster("processed_data/conservation_rasters/local_extirp.tif")
poss.decline <- raster("processed_data/conservation_rasters/poss_decline.tif")
pres.stable <- raster("processed_data/conservation_rasters/pres_stable.tif")



# Calculate null distribution of conservation status ----------------------
null_conserv_status <- function(pheno.raster) {
  data.frame(conservation = c("extirpated", "broad.extirp", "local.extirp", "poss.decline", "pres.stable"),
             area = c(sum(values(area(mask(extirpated, pheno.raster), na.rm = T)), na.rm = T),
                      sum(values(area(mask(broad.extirp, pheno.raster), na.rm = T)), na.rm = T),
                      sum(values(area(mask(local.extirp, pheno.raster), na.rm = T)), na.rm = T),
                      sum(values(area(mask(poss.decline, pheno.raster), na.rm = T)), na.rm = T),
                      sum(values(area(mask(pres.stable, pheno.raster), na.rm = T)), na.rm = T))) %>% 
    mutate(perc = area/sum(area))
}

null <- null_conserv_status(pheno.range)



# Calculate area of conservation status by color --------------------------
# Using the raster of predicted phenotypes (in prob(white)),
# calcualtes the area within each phenotypic class (specified by the thresholds)
# that falls within each area. 

# Optionally, also plots the distribution so you can look at it. 
conserv_status_by_color <- function(pheno.raster, white.thresh, brown.thresh) {
  
  # The phenotype rasters are given in prob(white), so will have to do 1-minus for prob of brown.
  
  # So, make three new rasters, one for each phenotype class, and then make them binary so
  # that the values there are either 1 or NA
  # "brown" area, where prob(brown) > brown thresh
  # "mixed" area, where white thresh < prob(brown) < brown thresh
  # "white" ares, where prob(brown) < white thresh
  brown <- pheno.raster
  values(brown) <- ifelse(1 - values(brown) > brown.thresh, 1, NA)
  
  white <- pheno.raster
  values(white) <- ifelse(1 - values(white) < white.thresh, 1, NA)
  
  # Mixed is just whatever isn't already in brown or white
  mixed <- mask(pheno.raster, brown, inverse = T) %>% 
    mask(., white, inverse = T)
  values(mixed) <- ifelse(!is.na(values(mixed)), 1, NA)
  
  
  # Set up the data frame
  # Order here must match that below
  phenos <- c("brown", "white", "mixed")
  cons.cats <- c("extirpated", "broad.extirp", "local.extirp", "poss.decline", "pres.stable") 
  results <- expand_grid(phenotype = phenos, conservation = cons.cats) %>% 
    mutate(area = NA)
  # A for loop to calculate it alldisk
  i <- 1
  for (pheno.cat in list(brown, white, mixed)) {
    for (cons.cat in list(extirpated, broad.extirp, local.extirp, poss.decline, pres.stable)) {
      results$area[i] <- sum(values(area(mask(pheno.cat, cons.cat), na.rm = T)), na.rm = T)
      i <- i + 1
    }
  }
  
  results %>% 
    pivot_wider(names_from = phenotype, values_from = area) 
}


# Calculate phenotypic status by area -------------------------------------
broad <- conserv_status_by_color(pheno.raster = pheno.range,
                                 white.thresh = 0.2,
                                 brown.thresh = 0.8)
narrow <- conserv_status_by_color(pheno.raster = pheno.range,
                                 white.thresh = 0.4,
                                 brown.thresh = 0.6)


# Do some stats -------------------------------------
broad.mat <- broad %>% 
  dplyr::select(-conservation) %>% 
  as.matrix(.)

broad.chisq.res <- tidy(chisq.test(broad.mat))

narrow.mat <- narrow %>% 
  dplyr::select(-conservation) %>% 
  as.matrix(.)

narrow.chisq.res <- tidy(chisq.test(narrow.mat))

# Calculating Cramer's V, an effect size stat for
# chi2 tests.
# Directons from https://www.ibm.com/support/knowledgecenter/SSEP7J_11.1.0/com.ibm.swg.ba.cognos.ug_ca_dshb.doc/cramersv.html

# Determine which field has the fewest number of categories. (color, with 3)
# Subtract 1 from the number of categories in this field. (So, 2)
# Multiply the result by the total number of records.
# Divide the chi-square value by the previous result. The chi-square value is obtained from the chi-square test of independence
# Take the square root.

total <- broad %>% 
  pivot_longer(brown:mixed, names_to = "phenos") %>% 
  summarize(total = sum(value))


cramerV.broad <- sqrt(broad.chisq.res$statistic/(total*2))
cramerV.narrow <- sqrt(narrow.chisq.res$statistic/(total*2))

# Save results -------------------------------------
save(null, broad, narrow, broad.chisq.res, narrow.chisq.res, cramerV.broad, cramerV.narrow, file = outfile)
