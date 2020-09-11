####################
# Project: WTJR
# Author: Timothy Thurman
# Purpose: Running SDMs with ENMeval
# Date Created: Fri Apr 24 12:39:58 2020
####################


# Load packages -----------------------------------------------------------
library(tidyverse)
# library(maptools)
# library(rgdal)
# library(maps)
library(raster)
library(dismo)
# library(maxnet)
library(ENMeval)
# library(stringr)
# # library(cowplot)
#library(rnaturalearth)
# source("src/misc_fxns.R")
# state_prov <- ne_states(c("united states of america", "canada"))

# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 csv of occurrence records for analysis
# 2 bioclim data
# 3 bg points
# 3 feature class
# 4 seed
# 5 cores
in_csv <- args[1]
in_bioclim <- args[2]
bg <- args[3]
fc <- args[4]
seed <- args[5]
cores <- args[6]

# in_csv <- "processed_data/thin/50km/wtjr_occ_50km_thin1.csv"
# in_bioclim <- "processed_data/bioclim_30arcsec_for_WTJR_SDM.tif"
# fc <- "L"
# seed <- 782
# cores <- 1
# bg <- "processed_data/bg_points_for_sdm.RData"
# Load occurence data -----------------------------------------------------
wtjr.occ <- read.csv(in_csv, stringsAsFactors = F) %>%
  dplyr::select(roundlon, roundlat)

# Load Bioclim data -------------------------------------------------------
bioclim <- brick(in_bioclim)
# Fix names
names(bioclim) <- paste("bio", 1:19, sep = "")
extent <- extent(bioclim)

# Background Points -----------------------------------------------------------------
load(bg)

# Parse input filename ----------------------------------------------------
dataset <- str_remove(in_csv, pattern = ".csv") %>% 
  str_split(., pattern = "_") %>% 
  .[[1]] %>%
  .[4:5] %>% 
  paste(., sep = "", collapse = "_")
  


# Loop to allow multiple running attempts ---------------------------------
set.seed(seed)
attempt_limit <- 1
target.file <- paste("results/enmeval_res_", dataset, "_", fc, ".RData", sep = "")


attempt <- 1
while (attempt <= (attempt_limit)) {
  # move the attempt counter up
  attempt <- attempt + 1
  if (!file.exists(target.file)) { # If the results file doesn't exist
    # Run the model within a try statement and save the results
    enmeval_res <- try(ENMevaluate(occ = wtjr.occ,
                                   env = bioclim,
                                   bg.coords = background,
                                   algorithm = "maxnet",
                                   method = "checkerboard2",
                                   parallel = F, rasterPreds = T, 
                                   numCores = cores,
                                   fc = fc))
    save(enmeval_res, file = target.file)
    
    # Then check to see if the result file is big enough
    if (file.size(target.file) < 1e2) { # if not
      # remove it and increment the attempt counter
      file.remove(target.file)
    } else {
      break
    }
  } else {
    break
  }
}

# # Collate results ---------------------------------------------------------
# # A loop to run on Griz to collect all the performance metrics
# # for each dataset and get them in a common data frame
# # From there, can choose which full ENMeval objects to download.
# 
# # Also, double-checks that the same background points were used in all analysis
# i <- 1
# for (dataset in c("wtjr_all_", "wtjr_thin1_", "wtjr_thin5_", "wtjr_thin10_", "wtjr_thin50_")) {
#   results <- NULL
#   for (file in list.files("results", full.names = T)) {
#     if (str_detect(string = file, pattern = dataset)) {
#       load(file)
#       # check background points the same as last one
#       if (i > 1) {
#         print(paste("For ", file, " do background points equal bg.points for last file? ", unique(last_one@bg.pts == enmeval_res@bg.pts), sep = ""))
#       }
#       results <- rbind(results, enmeval_res@results)
#       last_one <- enmeval_res
#       i <- i + 1
#     }
#   }
#   target.file <- paste("results/enmeval_", dataset, "metrics.RData", sep = "")
#   save(results, file = target.file)
# }
# 
# 
# 
# # Plot model evaluations --------------------------------------------------
# # Downloaded all results, then can do the following plotting locally:
# 
# 
# for (file in list.files("results/enmeval_round2_same_bg/", full.names = T)) {
#   load(file)
#   dataset <- str_split(file, pattern = "_")[[1]][6]
#   results %>% 
#     arrange(features, rm) %>% 
#     dplyr::select(features, rm, train.AUC, avg.test.AUC, avg.diff.AUC, avg.test.or10pct, AICc, parameters) %>% 
#     pivot_longer(train.AUC:parameters, names_to = "metric") %>% 
#     ggplot(aes(x = rm, y = value, group = features, color = features)) +
#     geom_point() +
#     geom_line() +
#     ggtitle(paste0("Performance metrics, dataset = ", dataset)) +
#     xlab("Regularization multiplier") +
#     ylab("performance metric value") +
#     facet_wrap(facets = vars(metric), scale = "free_y") +
#     theme_cowplot() +
#     ggsave(filename = paste0("results/enmeval_round2_same_bg/performance_plot_", dataset, ".pdf"), 
#            width = 12, height = 9, units = "in")
# }
# 
# 
# 
# 
# # Plot across datasets ----------------------------------------------------
# # Chose the best model from within each dataset (based on AICc), and can now compare across datasets
# # to see what to use as the final model
# for (file in list.files("results/enmeval_round2_same_bg/", full.names = T)) {
#   load(file)
#   dataset <- str_split(file, pattern = "_")[[1]][6]
#   results$dataset <- dataset
#   assign(x = paste0(dataset, "_metrics"), value = results)
#   rm(results)
# }
# 
# across_dataset_metrics <- rbind(all_metrics, thin1_metrics, thin5_metrics, thin10_metrics, thin50_metrics) %>% 
#   group_by(dataset) %>% 
#   filter(AICc == min(AICc)) %>% 
#   ungroup() %>% 
#   dplyr::select(dataset, features, rm, parameters, everything()) %>% 
#   dplyr::select(-contains("orMTP"),
#                  -contains("AIC"), -settings)
# 
# means <- across_dataset_metrics %>% 
#   dplyr::select(-contains("var")) %>% 
#   pivot_longer(cols = c(contains("avg"), train.AUC, parameters), names_to = "metric", values_to = "mean") %>% 
#   mutate(metric = str_remove(metric, pattern = "avg\\."))
# vars <- across_dataset_metrics %>% 
#   dplyr::select(-contains("avg"), -train.AUC, -parameters) %>% 
#   pivot_longer(cols = contains("var"), names_to = "metric", values_to = "vars") %>% 
#   mutate(metric = str_remove(metric, pattern = "var\\."))
# 
# 
# # Make and save plot
# means %>% 
#   left_join(vars) %>% 
#   mutate(sd = sqrt(vars),
#          n = 4) %>% 
#   mutate(se = sd/sqrt(n)) %>% 
#   mutate(metric = fct_relevel(metric, "train.AUC", "test.AUC",
#                               "diff.AUC", "test.or10pct", "parameters"),
#          dataset = fct_relevel(dataset, "all", "thin1", "thin5", "thin10", "thin50")) %>% 
#   ggplot(aes(x = dataset, y = mean, ymin = mean - 1.96*se, ymax = mean + 1.96*se, color = dataset)) +
#   geom_pointrange() +
#   ggtitle(paste0("Performance metrics across datasets")) +
#   xlab("dataset") +
#   ylab("value, +/- approx. 95% CI when possible") +
#   facet_wrap(facets = vars(metric), scale = "free_y") +
#   theme_cowplot() +
#   ggsave(filename = "results/enmeval_round2_same_bg/performance_across_datasets.pdf", 
#          width = 12, height = 7, units = "in")
# 
# means %>% 
#   left_join(vars) %>% 
#   mutate(sd = sqrt(vars),
#          n = 4) %>% 
#   mutate(se = sd/sqrt(n)) %>% 
#   mutate(metric = fct_relevel(metric, "train.AUC", "test.AUC",
#                               "diff.AUC", "test.or10pct", "parameters"),
#          dataset = fct_relevel(dataset, "all", "thin1", "thin5", "thin10", "thin50")) %>% 
#   ggplot(aes(x = dataset, y = mean, ymin = mean - 1.96*se, ymax = mean + 1.96*se, color = dataset)) +
#   geom_pointrange() +
#   ggtitle(paste0("Performance metrics across datasets")) +
#   xlab("dataset") +
#   ylab("value, +/- approx. 95% CI when possible") +
#   facet_wrap(facets = vars(metric), scale = "free_y") +
#   theme_cowplot() +
#   ggsave(filename = "results/enmeval_round2_same_bg/performance_across_datasets.png", 
#          width = 12, height = 7, units = "in")
#   
# 
# # From these data, the all data dataset seems best.
# # Let's download it and take a look at the output. 
# 
# 
# # SDM from all data -------------------------------------------------------
# load("results/enmeval_round2_same_bg/enmeval_res_wtjr_all_LQHPT_check2.RData")
# 
# # Extract predictions from best model
# best.preds <- subset(enmeval_res@predictions, "LQHPT_1")
# 
# # evaluate the best model
# eval.best.mod <- dismo::evaluate(p = raster::extract(x = best.preds, y = enmeval_res@occ.pts), 
#                                  a = raster::extract(x = best.preds, y = enmeval_res@bg.pts))
# # Use evaluation to get a threshold
# spec_sens_thresh <- dismo::threshold(eval.best.mod, stat = "spec_sens")
# # Get sensitivity at spec_sens thresh
# eval.best.mod@TPR[eval.best.mod@t == spec_sens_thresh]
# # Get specificity at spec_sens thresh
# eval.best.mod@TNR[eval.best.mod@t == spec_sens_thresh]
# 
# eval.best.mod@TNR[eval.best.mod@t == sens_95_thresh]
# 
# 
# 
# sens_95_thresh <- dismo::threshold(eval.best.mod, stat = "sensitivity", sensitivity = 0.95)
# sens_99_thresh <- dismo::threshold(eval.best.mod, stat = "sensitivity", sensitivity = 0.99)
# 
# 
# # Make the range raster and plot it, saving as both a png and a pdf
# range.raster.spec_sens <- best.preds > spec_sens_thresh
# range.raster.sens95 <- best.preds > sens_95_thresh
# range.raster.sens99 <- best.preds > sens_99_thresh
# 
# # Save as a binary pres/absence file --------------------------------------
# values(range.raster.spec_sens) <- ifelse(values(range.raster.spec_sens) < 0.5, NA, 1)
# values(range.raster.sens95) <- ifelse(values(range.raster.sens95) < 0.5, NA, 1)
# values(range.raster.sens99) <- ifelse(values(range.raster.sens99) < 0.5, NA, 1)
# 
# writeRaster(range.raster.spec_sens, filename = "results/enmeval_round2_same_bg/range_alldata_best_specsens", overwrite = T)
# writeRaster(range.raster.sens95, filename = "results/enmeval_round2_same_bg/range_alldata_best_sens95", overwrite = T)
# # writeRaster(range.raster.sens99, filename = "results/enmeval_round2_same_bg/range_alldata_best_specsens", overwrite = T)
