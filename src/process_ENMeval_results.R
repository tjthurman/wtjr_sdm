####################
# Project: WTJR
# Author: Timothy Thurman
# Purpose: Processing ENMeval results
# Date Created: Thu Spet 24 11:22:58 2020
####################


# Load packages -----------------------------------------------------------
library(tidyverse)
library(raster)
library(dismo)
library(ENMeval)


# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# # 1 folder of enmeval results
# in_csv <- args[1]
# in_bioclim <- args[2]
# bg <- args[3]
# fc <- args[4]
# seed <- args[5]
# cores <- args[6]

res.folder <- args[1]


# Load results and plot them ----------------------------------------------
# thin_dists <- paste0(c(0,1, 5,10,50), "km")

results <- NULL
for (file in list.files(res.folder, full.names = T)) {
  # Load all enmeval results objects, collate their metrics into one data frame
  if (str_detect(string = file, pattern = "enmeval_res_")) {
    dataset <- str_extract(string = file, pattern = "\\d+km")
    load(file)
    oneres <- enmeval_res@results
    oneres$thin_dist <- dataset
    results <- rbind(results, oneres)
    rm(enmeval_res, oneres)
  }
}
# Write those results out to a file
sorted <- results %>% 
  dplyr::select(thin_dist, everything()) %>%
  dplyr::select(-settings) %>% 
  arrange(match(thin_dist, c("0km", "1km", "5km", "10km", "50km")), 
          match(features, c("L", "LQ", "H", "LQH", "LQHP", "LQHPT")), rm) 

write.csv(x = sorted, file = "results/enmeval/enmeval_metrics.csv", row.names = F)

# Make a plot for each dataset
for (dist in unique(sorted$thin_dist)) {
  sorted %>%
    filter(thin_dist == dist) %>% 
    dplyr::select(features, rm, train.AUC, avg.test.AUC, avg.diff.AUC, avg.test.or10pct, AICc, w.AIC, parameters) %>%
    pivot_longer(train.AUC:parameters, names_to = "metric") %>%
    ggplot(aes(x = rm, y = value, group = features, color = features)) +
    geom_point() +
    geom_line() +
    ggtitle(paste0("Performance metrics, dataset = ", dist)) +
    xlab("Regularization multiplier") +
    ylab("performance metric value") +
    facet_wrap(facets = vars(metric), scale = "free_y") +
    ggsave(filename = paste0("results/enmeval/performance_plot_", dist, ".pdf"),
           width = 12, height = 9, units = "in")
}
  
# Get .csv of the best parameters for each dataset. 
best_mods <- sorted %>% 
  group_by(thin_dist) %>% 
  filter(AICc == min(AICc))
write.csv(x = best_mods, file = "results/enmeval/enmeval_best_model_per_thin_AIC.csv", row.names = F)


# Performance plot of best model in each dataset
means <- best_mods %>%
  dplyr::select(-contains("var")) %>%
  pivot_longer(cols = c(contains("avg"), train.AUC, parameters), names_to = "metric", values_to = "mean") %>%
  mutate(metric = str_remove(metric, pattern = "avg\\."))
vars <- best_mods %>%
  dplyr::select(-contains("avg"), -train.AUC, -parameters) %>%
  pivot_longer(cols = contains("var"), names_to = "metric", values_to = "vars") %>%
  mutate(metric = str_remove(metric, pattern = "var\\."))
 
means %>%
  left_join(vars) %>%
  mutate(sd = sqrt(vars),
         n = 4) %>%
  mutate(se = sd/sqrt(n)) %>%
  mutate(metric = fct_relevel(metric, "train.AUC", "test.AUC",
                              "diff.AUC", "test.or10pct", "parameters")) %>% 
  separate(thin_dist, into = c("thin_dist", "extra"), sep = -2, convert = T) %>%
  ggplot(aes(x = as.factor(thin_dist), y = mean, ymin = mean - 1.96*se, ymax = mean + 1.96*se, color = as.factor(thin_dist))) +
  geom_pointrange() +
  ggtitle(paste0("Performance metrics across datasets")) +
  xlab("Thinning distance of dataset") +
  ylab("value, +/- approx. 95% CI when possible") +
  facet_wrap(facets = vars(metric), scale = "free_y") +
  ggsave(filename = "results/enmeval/performance_plot_best_models.pdf",
         width = 12, height = 7, units = "in")
