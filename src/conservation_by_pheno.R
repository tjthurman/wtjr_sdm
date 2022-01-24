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
library(cowplot)

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


# Proportion of mismatch by conservation status ---------------------------
# Jeff had the idea to more-or-less repeat the analysis we do above, 
# comparing conservation status and current phenotype,
# but instead to compare conservation status and predicted mismatch (map from fig. 4)


# So, first we can just calculate mismatch, copying the code that we use to do it in Figure 4:

# Get the files that contain the current and future predicted probWhite based on SRT
current_file <- "results/pheno/current_predicted_probWhite_SDMrange_SRT.tif"
future_file <- "results/pheno/future_predicted_probWhite_SDMrange.tif"
future_pheno <- raster(future_file)
current_pheno <- raster(current_file)

# Then calculate change in probBrown from those maps
mismatch <- (1- future_pheno) - (1 - current_pheno) 


# Now, the area where we can calculate mismatch is smaller than the full range
# we can use for the current phenotypes (doesn't include most of Canada)

# So, first, we have to re-calculate the distribution of how much
# area is in each conservation category. To be honest, now that I do this I'm not sure we need
# this calculation (or the calculation above for the full range): I don't think we use it again??
# but will keep for now, and edit out if we don't need it anymore. 


# Can mostly just reuse the code in the function above, with an extra crop
null_SRT_range <- data.frame(conservation = c("extirpated", "broad.extirp", "local.extirp", "poss.decline", "pres.stable"),
                             area = c(sum(values(area(mask(crop(extirpated, mismatch), mismatch), na.rm = T)), na.rm = T),
                                      sum(values(area(mask(crop(broad.extirp, mismatch), mismatch), na.rm = T)), na.rm = T),
                                      sum(values(area(mask(crop(local.extirp, mismatch), mismatch), na.rm = T)), na.rm = T),
                                      sum(values(area(mask(crop(poss.decline, mismatch), mismatch), na.rm = T)), na.rm = T),
                                      sum(values(area(mask(crop(pres.stable, mismatch), mismatch), na.rm = T)), na.rm = T))) %>% 
  mutate(perc = area/sum(area))

# Just as an eyeball test, compare the two
null
null_SRT_range
# The SRT range (which doesn't include Canada), has basically no areas that are presumed stable, which makes sense:
# that is only Manitoba. 


# Now, on to some further analysis
# First, we can just make a violin plot: x-axis as the conservation categories, y axis as mismatch

# To do that, I need to extract all the mismatch values for a given area for each conservation status
# Not the nicest way to do this, I suppose, but will just loop over a vector of the names of the objects
# for each conservation status.
# Then, for each conservation status, mask the mismatch against the conservation status,
# with an extra crop to get the extents to map. 

# This will output one large data frame that I can use for plotting the mismatch values by category. 
mismatch_values_by_cons <- tibble()
for (cons.cat in c("extirpated", "broad.extirp", "local.extirp", "poss.decline", "pres.stable")) {
  single <- tibble(conserv_status = cons.cat,
                   mismatch = values(mask(mismatch, crop(get(cons.cat), mismatch)))) %>% 
    filter(!is.na(mismatch))
  
  mismatch_values_by_cons <- rbind(mismatch_values_by_cons, single)                                    
}

# Do a little re-factoring to make the plot a little more understandable
mismatch_values_by_cons$conserv_status <- factor(mismatch_values_by_cons$conserv_status, 
                             levels = c("extirpated", "broad.extirp", "local.extirp", "poss.decline", "pres.stable"),
                             labels = c("extirpated", "broad\nextirpations", "local\nextirpations",
                                        "possible\ndeclines", "presumed\nstable"))
# And then plot
mismatch_values_by_cons %>%
  ggplot(aes(x = conserv_status, y = mismatch, fill = conserv_status)) +
  geom_violin(scale = "width", draw_quantiles = c(0.5)) +
  scale_fill_manual(
    values = c(
      "extirpated" = "#EC4741",
      "broad\nextirpations" = "#F8B732",
      "local\nextirpations" = "#F2F031",
      "possible\ndeclines" = "#ABD344",
      "presumed\nstable" = "#44854B"
    ),
    name = "conservation\nstatus",
  ) +
  guides(fill = "none") +
  theme(legend.position = "none") +
  theme_cowplot()

# So, just by eyeballing we can see that future mismatch is non-randomly distributed across the conservation categories
# Interestingly, this is perhaps more related to current color than mismatch itself.
# That is, the two conservation categories with the lowest average mismatch are "extirpated" and "presumed stable",
# but for different reasons.

# Extirpated areas are are mostly already brown and in the south, so they have little future mismatch.
# presumed stable is really only a very small part of southern Manitoba that is included in the SRT dataset
# And it is aleady white and likely to stay whie, to no mismatch. 

# Within the other conservation categories, the distributions are actually relatively similar, if generally bimodal:
# some part that is low mismatch and some part that is high mismatch. But, I expect that the regions that are low mismatch
# have different causes across categories: possible declines are more in the north, those low mismatch areas are low mismatch because
# they're white and will stay white. Broad extirpations are more in the south, and their low mismatch are low mismatch because
# they are already brown and will stay brown. 

# In any case, we can try to put some stats on this variation.
# One option might be to do a bunch of pairwise KS tests: is there a procedure for that?

# In any case, an easier option would be to just do the chisq test, 
# as we did with phenotype by conservation status.
# to do that, we need to turn our continuous
# mismatch measure into discretized categories. 

# We can of course choose different bins to do this. 
# Here are a couple options: one with 3 categories, and one with 4
mismatch_disc_3 <- cut(mismatch, breaks = c(0, 0.3, 0.6, 1))
mismatch_disc_4 <- cut(mismatch, breaks = c(0, 0.1, 0.3, 0.5, 1))
plot(mismatch_disc_3)
plot(mismatch_disc_4)

# Just eyeballing, I think 4 categories works nicer.
# We could think of it being these 4 categories: 

# "no" mismatch is 0-0.1
# "minor" mismatch is 0.1-0.3
# "moderate" mismatch is 0.3-0.5
# "extreme" mismatch is 0.5 to 1 (really, to 0.85, the max we see)

# To do the chi-sq test, we now need to do the same thing that we did above with phenotypes:
# loop over both the mismatch categories and the conservation categories,
# and finding the area of the overlapping rasters for each of them.
# That is, finding the area of extreme mismatch that is in extirpated, and in broad extirpations, etc. 

# Again, will modify the loop, but will do it maybe a little nicer, looping over the 
# names of the categories. 
# The raster cut I did above is nice, in that it assigns integers to each of the categories,
# so I'll take advantage of that for subsetting. 

# The below code only works for this 4-category classification of mismatch,
# would need to do some minor edits to make it work with 3 categories. 
mismatch_cat <- c("none", "minor", "moderate", "extreme")
cons.cats <- c("extirpated", "broad.extirp", "local.extirp", "poss.decline", "pres.stable") 
results <- expand_grid(mismatch_cat = mismatch_cat, conservation = cons.cats) %>% 
  mutate(area = NA,
         check = NA)


i <- 1
for (mismatch_int in 1:length(mismatch_cat)) {
  # Take the discretized mismatch from above
  mismatch_raster_disc <- mismatch_disc_4
  # Redo-values to make it one group at a time: 1 if it belongs to the group,
  # NA if it doesn't.
  # Equivalent to masking out the other categories
  values(mismatch_raster_disc) <- ifelse(values(mismatch_raster_disc) == mismatch_int, 1, NA)
  # Then, for each of the above mismatch categories, loop over the conservation categories
  for (cons in cons.cats) {
    # And find the area where the two categories overlap
    results$area[i] <- sum(values(area(mask(mismatch_raster_disc, crop(get(cons), mismatch_raster_disc)), na.rm = T)), na.rm = T)
    results$check[i] <- paste(mismatch_cat[mismatch_int], cons)
    i <- i + 1
  }
}

# Then get it into the matrix format we need for the chi-sq test
mismatch_disc_area_mat <- results %>% 
  dplyr::select(-check) %>% 
  pivot_wider(names_from = mismatch_cat, values_from = area) %>% 
  dplyr::select(-conservation) %>% 
  as.matrix(.)

# And do the chi-sq test
mismatch_disc.chisq.res <- tidy(chisq.test(mismatch_disc_area_mat))
# Yes, non-randomly distributed
# But, this is maybe to be expected: the "N" for this is gigantic,
# as it is the total area of the range of WTJR in the US (the statistical unit is 1 sq km)
# So, even biologically insiginificant variation is likely to be statistically significant. 

# So, as above, we can use Cramer's V to get an effect size 
total_mismatch <- results %>% 
  summarize(total = sum(area))
cramerV.mismatch_disc <- sqrt(mismatch_disc.chisq.res$statistic/(total_mismatch*3))
# So, a smaller effect than the result for current color

# Wikipedia mentions that Cramer's V can be a biased estimator, and suggests using a bias correction
# Looks like there's also an R package that supplies a cramers V calculation.
# So I can try that:
# renv::install("rcompanion")
library(rcompanion)
cramerV(mismatch_disc_area_mat, bias.correct = F) # .1241, same as what I calculated
cramerV(mismatch_disc_area_mat, bias.correct = T) # Bias correction made no difference

# Can double-check the others, just to be sure:
cramerV(broad.mat, bias.correct = F)
cramerV(broad.mat, bias.correct = T)

cramerV(narrow.mat, bias.correct = F)
cramerV(narrow.mat, bias.correct = T)
# So I get the same results my way as with rcompanion, and the bias correction makes no difference. 








