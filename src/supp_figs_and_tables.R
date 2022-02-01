####################
# Project: WTJR SDMs
# Author: Timothy Thurman
# Purpose: Figures and tables for supplemental material
# Date Created: Fri Oct 21 10:23:30 2020
####################


# Load packages -----------------------------------------------------------
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(raster)
library(ggsn)
library(broom)
library(gridExtra)
library(cowplot)
library(lemon)
library(RColorBrewer)



# # for running in a pipline
args = commandArgs(trailingOnly=TRUE)
# 1 snow cover glm rdata
# 2 srt glm rdata
# 3 conservation pheno overlap
# 4 current SRT predicted phenos
# 5 future SRT predicted phenos
# 6 current snow cover predicted phenos
# 7 csv of compiled sims, additive 2 locus consK
# 8 csv of compiled sims, recessive 2 locus consK
# 9 csv of compiled sims, additive 2 locus varyK
# 10 csv of compiled sims, recessive 1 locus consK

# 11 snow cover table out
# 12 snow cover metrics out
# 13 srt table out
# 14 srt metrics out
# 15 broad chisq out
# 16 pheno comparison map out
# 17 glm method compare map out
# 18 percent brown change map out
# 19 discrete current pheno map out
# 20 extended data arch fig out, pdf
# 21 extended data arch fig out, jpeg
# 22 extended data robust fig out, pdf
# 24 extended data robust fig out, jpeg


# Inputs
snow_cover_glm <- args[1]
srt_glm <- args[2]
conservation_overlap <- args[3]
current_file <- args[4]
future_file <- args[5]
current_cover_file <- args[6]
best_metrics_file <- args[7]
additive_consK_2locus_late_file <- args[8]
recessive_consK_2locus_late_file <- args[9]
additive_varyK_2locus_late_file <- args[10]
recessive_consK_1locus_late_file <- args[11]

# Outputs
snow_cover_table <- args[12]
snow_cover_metrics <- args[13]
srt_table <- args[14]
srt_metrics  <- args[15]
broad_chisq_out <- args[16]
pheno_compare_maps <- args[17]
ext_data_sdm_fig <- args[18]
ext_data_sdm_fig_jpg <- args[19]
glm_method_compare_metrics <- args[20]
arch_varyK_fig_pdf <- args[21]
arch_varyK_fig_jpeg <- args[22]
robust_fig_pdf <- args[23]
robust_fig_jpeg <- args[24]


# For running as a script
# Inputs
# snow_cover_glm <- "results/pheno/current_pheno_glm.RData"
# srt_glm <- "results/pheno/current_pheno_glm_SRT.RData"
# conservation_overlap <- "results/conservation/cons_by_current_color.RData"
# current_file <- "results/pheno/current_predicted_probWhite_SDMrange_SRT.tif"
# future_file <- "results/pheno/future_predicted_probWhite_SDMrange.tif"
# current_cover_file <- "results/pheno/current_predicted_probWhite_SDMrange.tif"
# best_metrics_file <- "results/enmeval/enmeval_best_model_per_thin_AIC.csv"
# additive_consK_2locus_late_file <- "results/slim_summaries/additive_constantK_2locus_late.csv"
# recessive_consK_2locus_late_file <- "results/slim_summaries/recessive_constantK_2locus_late.csv"
# additive_varyK_2locus_late_file <- "results/slim_summaries/additive_varyK_2locus_late.csv"
# recessive_consK_1locus_late_file <- "results/slim_summaries/recessive_constantK_1locus_late.csv"
# 
# 
# # outputs
# snow_cover_table <- "results/pheno/glm_table_current_snow_cover.csv"
# snow_cover_metrics <- "results/pheno/glm_metrics_current_snow_cover.csv"
# srt_table <- "results/pheno/glm_table_current_srt.csv"
# srt_metrics <- "results/pheno/glm_metrics_current_srt.csv"
# broad_chisq_out <- "results/conservation/broad_chisq_res.csv"
# ext_data_sdm_fig <- "results/figures/supplemental/extended_data_SDMs.pdf"
# ext_data_sdm_fig_jpg <- "results/figures/supplemental/extended_data_SDMs.jpeg"
# glm_method_compare_metrics <- "results/pheno/model_difference_metrics.csv"
# arch_varyK_fig_pdf <- "results/figures/supplemental/extended_data_arch_varyK_sim_res.pdf"
# arch_varyK_fig_jpeg <- "results/figures/supplemental/extended_data_arch_varyK_sim_res.jpeg"
# robust_fig_pdf <- "results/figures/supplemental/extended_data_sim_robust.pdf"
# robust_fig_jpeg <- "results/figures/supplemental/extended_data_sim_robust.jpeg"


# Load map data -----------------------------------------------------------
state_prov <- ne_states(c("united states of america", "canada"), returnclass = "sf")
countries <- ne_countries(continent = "north america", scale = 10, returnclass = "sf")
coast <- ne_coastline(scale = 10, returnclass = "sf")


# TABLES ------------------------------------------------------------------
# GLM tables --------------------------------------------------------------
load(snow_cover_glm)
load(srt_glm)

tidy(lepto.mod) %>% 
  write.csv(., snow_cover_table, row.names = F)
glance(lepto.mod) %>% 
  mutate(pseudo.r2 = 1 - (deviance/null.deviance)) %>% 
  write.csv(., snow_cover_metrics, row.names = F)

tidy(srt.mod) %>% 
  write.csv(., srt_table, row.names = F)
glance(srt.mod) %>% 
  mutate(pseudo.r2 = 1 - (deviance/null.deviance)) %>% 
  write.csv(., srt_metrics, row.names = F)


# Chi2 and cramers v results ----------------------------------------------
load(conservation_overlap)

broad.chisq.res %>% 
  mutate(cramerV = cramerV.broad[[1]]) %>% 
  write.csv(., broad_chisq_out, row.names = F)


# FIGURES -----------------------------------------------------------------
# Load predicted phenotypes ----------------------------------------------------
future_pheno <- raster(future_file) %>% 
  aggregate(., fact = 4, expand = F)
current_pheno <- raster(current_file) %>% 
  aggregate(., fact = 4, expand = F)

# Plots of current and future prob(brown) ---------------------------------
current_probBrown_df <- as(current_pheno, "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1])) %>% 
  mutate(probBrown = 1 - current_predicted_probWhite_SDMrange_SRT)
names(current_probBrown_df) <- c("probWhite", "Long", "Lat", "probBrown")

future_probBrown_df <- as(future_pheno, "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1])) %>% 
  mutate(probBrown = 1 - future_predicted_probWhite_SDMrange)
names(future_probBrown_df) <- c("probWhite", "Long", "Lat", "probBrown")

# Color pallette
pal <- colorRampPalette(colors = rev(c("#69431A", "floralwhite"))) # Colors from Mafalda's figs

map_pheno_current <- ggplot() +
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255), size = 0.25) +
  geom_tile(data = current_probBrown_df, aes(fill = probBrown, color = probBrown, x = Long, y = Lat)) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, size = 0.25) +
  geom_sf(data = countries, color = "grey10", fill = NA, size = 0.25) +
  geom_sf(data= coast, color = "grey10", fill = NA, size = 0.25) +
  coord_sf(
    xlim = c(-125, -89.5),
    ylim = c(32, 49.5),
    clip = "on", 
    expand = F) +
  scale_fill_gradientn(colors = pal(100), guide = F) +
  scale_color_gradientn(colors = pal(100), name = NULL) +
  guides(color = guide_colorbar(label.position = "right", ticks = T,
                                ticks.colour = "black", frame.colour = "black", 
                                frame.linewidth = 0.75, 
                                barwidth = unit(1.75, "mm"), 
                                barheight = unit(8, "mm"))) +
  scale_y_continuous(breaks = seq(34, 48, by = 4)) + 
  scale_x_continuous(breaks = seq(-120, -90, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.88, 0.24),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.ticks = element_line(size = 0.75),
        axis.text = element_text(size = 5, color = "black"),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 0.75),
        plot.margin = unit(c(1, 1, 1, 1), "mm"),
        plot.title = element_text(size = 7, margin = margin(0,0,0,0, "mm"))) +
  scalebar(x.min = -125, x.max = -89.5,
           y.min = 32, y.max = 49.5, 
           dist  = 250, dist_unit = "km", model = "WGS84", transform = T, 
           anchor = c(x = -92.15, y = 33.3),
           st.size = 5/.pt, border.size = 0.25) +
  ggtitle("Current Prob(Brown), SRT model")


map_pheno_future <- ggplot() +
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255), size = 0.25) +
  geom_tile(data = future_probBrown_df, aes(fill = probBrown, color = probBrown, x = Long, y = Lat)) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, size = 0.25) +
  geom_sf(data = countries, color = "grey10", fill = NA, size = 0.25) +
  geom_sf(data= coast, color = "grey10", fill = NA, size = 0.25) +
  coord_sf(
    xlim = c(-125, -89.5),
    ylim = c(32, 49.5),
    clip = "on", 
    expand = F) +
  scale_fill_gradientn(colors = pal(100), guide = F) +
  scale_color_gradientn(colors = pal(100), name = NULL) +
  guides(color = guide_colorbar(label.position = "right", ticks = T,
                                ticks.colour = "black", frame.colour = "black", 
                                frame.linewidth = 0.75, 
                                barwidth = unit(1.75, "mm"), 
                                barheight = unit(8, "mm"))) +
  scale_y_continuous(breaks = seq(34, 48, by = 4)) + 
  scale_x_continuous(breaks = seq(-120, -90, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.88, 0.24),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.ticks = element_line(size = 0.75),
        axis.text = element_text(size = 5, color = "black"),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 0.75),
        plot.margin = unit(c(1, 3, 1, 1), "mm"),
        plot.title = element_text(size = 7, margin = margin(0,0,0,0, "mm"))) +
  scalebar(x.min = -125, x.max = -89.5,
           y.min = 32, y.max = 49.5, 
           dist  = 250, dist_unit = "km", model = "WGS84", transform = T, 
           anchor = c(x = -92.15, y = 33.3),
           st.size = 5/.pt, border.size = 0.25) +
  ggtitle("Future Prob(Brown), SRT model")



# Difference in SRT vs. snow cover models ---------------------------------
# first, calculate model difference metrics on the full
# resolution data

# First, get the current SRT results, but do not aggregate
cur_SRT_full <- raster(current_file) 
# Then, get current snow cover without aggregation
cur_SC_full <- raster(current_cover_file)

# reproject the current SRT to be in terms of the current SC
cur_SRT_full <- projectRaster(from = cur_SRT_full, to = cur_SC_full)

# mask the current snow cover
cur_SC_full_masked <- mask(cur_SC_full, cur_SRT_full)


# Get model difference, do it in terms of probBrown
model_difference <- (1- cur_SC_full_masked) - (1 - cur_SRT_full)

# Convert to DF
model_difference_current_probBrown <- model_difference %>% 
  as(., "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1])) 
names(model_difference_current_probBrown) <- c("difference", "Long", "Lat")


# Calculate and save some metrics
model_diff_metrics <- tibble(mean_diff = mean(model_difference_current_probBrown$difference),
                             median_diff = median(model_difference_current_probBrown$difference),
                             sd_diff = sd(model_difference_current_probBrown$difference),
                             MAE_diff = mean(abs(model_difference_current_probBrown$difference))) %>% 
  write.csv(., glm_method_compare_metrics, row.names = F)

# Then, aggregate to a lower resolution for plotting
model_difference_lowres <- aggregate(model_difference, fact = 4)

# convert to DF
model_difference_current_probBrown_lowres <- model_difference_lowres %>% 
  as(., "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1])) 
names(model_difference_current_probBrown_lowres) <- c("difference", "Long", "Lat")

# Plot
diff_pal <- colorRampPalette(colors = c("#762a83", "#f7f7f7", "#1b7837"))

map_model_difference <- ggplot() +
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255), size = 0.25) +
  geom_tile(data = model_difference_current_probBrown_lowres, aes(fill = difference, color = difference, x = Long, y = Lat)) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, size = 0.25) +
  geom_sf(data = countries, color = "grey10", fill = NA, size = 0.25) +
  geom_sf(data= coast, color = "grey10", fill = NA, size = 0.25) +
  coord_sf(
    xlim = c(-125, -89.5),
    ylim = c(32, 49.5),
    clip = "on", 
    expand = F) +
  scale_fill_gradientn(colors = diff_pal(100), guide = F) +
  scale_color_gradientn(colors = diff_pal(100), name = NULL) +
  guides(color = guide_colorbar(label.position = "right", ticks = T,
                                ticks.colour = "black", frame.colour = "black", 
                                frame.linewidth = 0.75, 
                                barwidth = unit(1.75, "mm"), 
                                barheight = unit(8, "mm"))) +
  scale_y_continuous(breaks = seq(34, 48, by = 4)) + 
  scale_x_continuous(breaks = seq(-120, -90, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.88, 0.24),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.ticks = element_line(size = 0.75),
        axis.text = element_text(size = 5, color = "black"),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 0.75),
        plot.margin = unit(c(1, 3, 1, 1), "mm"),
        plot.title = element_text(size = 7, margin = margin(0,0,0,0, "mm"))) +
  scalebar(x.min = -125, x.max = -89.5,
           y.min = 32, y.max = 49.5, 
           dist  = 250, dist_unit = "km", model = "WGS84", transform = T, 
           anchor = c(x = -92.15, y = 33.3),
           st.size = 5/.pt, border.size = 0.25) +
  ggtitle("Difference in current Prob(Brown) between snow cover and SRT models")



# Map of current, discrete phenos -----------------------------------------
pred.new <- raster(current_cover_file) %>% 
  aggregate(., fact = 4, expand = F) 

pred.pheno.df <- as(pred.new, "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1]))
names(pred.pheno.df) <- c("probWhite", "Long", "Lat")

pred.pheno.df$discretePheno <- NA
pred.pheno.df$discretePheno[pred.pheno.df$probWhite < 0.2] <- "brown"
pred.pheno.df$discretePheno[pred.pheno.df$probWhite > 0.8] <- "white"
pred.pheno.df$discretePheno[pred.pheno.df$probWhite < 0.8 & pred.pheno.df$probWhite > 0.2] <- "mixed"

discrete_pheno <- ggplot() +
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255)) +
  geom_tile(data = pred.pheno.df, aes(fill = discretePheno, color = discretePheno, x = Long, y = Lat)) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, ) +
  geom_sf(data = countries, color = "grey10", fill = NA) +
  geom_sf(data= coast, color = "grey10", fill = NA) +
  coord_sf(
    xlim = c(-132, -89),
    ylim = c(32, 56.5),
    clip = "on", 
    expand = F) +
  scale_fill_manual(values = c(
    "brown" = "#69431A",
    "mixed" = "#9F8567",
    "white" = "floralwhite"
  ), name = "phenotype") +
  scale_color_manual(values = c(
    "brown" = "#69431A",
    "mixed" = "#B39D83", # midpoint of our palette
    "white" = "floralwhite"
  ), name = "phenotype", guide = F) +
  scale_y_continuous(breaks = seq(35, 55, by = 5)) + 
  scale_x_continuous(breaks = seq(-130, -70, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.88, 0.15),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.ticks = element_line(size = 0.75),
        axis.text = element_text(size = 5, color = "black"),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 0.75),
        plot.margin = unit(c(1, 1, 1, 1), "mm"),
        plot.title = element_text(size = 7, margin = margin(0,0,0,0, "mm"))) +
  ggsn::scalebar(x.min = -132, x.max = -80.5,
                 y.min = 32, y.max = 56.5, 
                 dist  = 300, dist_unit = "km", model = "WGS84", transform = T, 
                 anchor = c(x = -122.48, y = 35),
                 st.size = 5/.pt, border.size = 0.25) +
  ggtitle("Current Prob(Brown), snow cover model, discretized")



# SDM performance metrics -------------------------------------------------
best_mods <- read.csv(best_metrics_file)

# Performance plot of best model in each dataset
means <- best_mods %>%
  dplyr::select(-contains("var")) %>%
  pivot_longer(cols = c(contains("avg"), train.AUC, parameters), names_to = "metric", values_to = "mean") %>%
  mutate(metric = str_remove(metric, pattern = "avg\\."))
vars <- best_mods %>%
  dplyr::select(-contains("avg"), -train.AUC, -parameters) %>%
  pivot_longer(cols = contains("var"), names_to = "metric", values_to = "vars") %>%
  mutate(metric = str_remove(metric, pattern = "var\\."))

sdm_metrics_plot <- means %>%
  left_join(vars) %>%
  mutate(sd = sqrt(vars),
         n = 4) %>%
  mutate(se = sd/sqrt(n)) %>%
  mutate(metric = fct_relevel(metric, "train.AUC", "test.AUC",
                              "diff.AUC", "test.or10pct", "parameters")) %>% 
  separate(thin_dist, into = c("thin_dist", "extra"), sep = -2, convert = T) %>%
  ggplot(aes(x = as.factor(thin_dist), y = mean, ymin = mean - 1.96*se, ymax = mean + 1.96*se)) +
  geom_pointrange(size = 0.1) +
  ggtitle(paste0("Performance metrics across datasets")) +
  xlab("Thinning distance of dataset (km)") +
  ylab("value, +/- approx. 95% CI when possible") +
  facet_rep_wrap(facets = vars(metric), scales = "free_y", nrow = 2, ncol = 3) +
  theme_cowplot() +
  theme(axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        plot.title =  element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        strip.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")))


# Percent of range that is brown ------------------------------------------
len.vec.current <- length(current_probBrown_df$probBrown)
len.vec.future <- length(future_probBrown_df$probBrown)


thresholds <- seq(.2, .8, length.out = 200)
area_compare <- data.frame(prob.brown.thresh = thresholds,
                           current = NA,
                           future =NA )
i <- 1
for (br.thresh in thresholds) {
  area_compare$current[i] <- sum(current_probBrown_df$probBrown > br.thresh, na.rm = T)/len.vec.current
  area_compare$future[i] <- sum(future_probBrown_df$probBrown > br.thresh, na.rm = T)/len.vec.future
  i <- i + 1
}

perc_brown_plot <- area_compare %>% 
  pivot_longer(current:future, names_to = "time_period", values_to = "perc.brown") %>% 
  ggplot(aes(x = prob.brown.thresh, y = perc.brown, color = time_period)) +
  geom_line() +
  xlab("Threshold of Prob(brown)") +
  ylab("Percent of range > threshold (\"brown\")") +
  scale_color_discrete(type = c(alpha("purple", 1), alpha("darkgreen", 1))) +
  theme_cowplot() + 
  theme(legend.position = c(0.75, 0.85),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.title = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(2, 3, 2, 2), "mm")) 



# Combine SDM figs into one large figure ----------------------------------
top_row <- cowplot::plot_grid(sdm_metrics_plot, perc_brown_plot,
                              nrow = 1, ncol = 2,
                              labels = c("a", "b"),
                              rel_widths = c(2, 1),
                              label_fontfamily = "sans", 
                              label_fontface = "bold", 
                              label_size = 8,
                              align = "none")

bottom_rows <- cowplot::plot_grid(discrete_pheno, map_model_difference,
                            map_pheno_current, map_pheno_future,
                            nrow = 2, ncol = 2, 
                            labels = c("c", "d", "e", "f"), 
                            label_fontfamily = "sans", 
                            label_fontface = "bold", 
                            label_size = 8,
                            align = "none")

combo <- cowplot::plot_grid(top_row, bottom_rows,
                           nrow = 2, ncol = 1,
                           labels = NA,
                           rel_heights = c(0.3, 0.7))


ggsave(plot = combo, ext_data_sdm_fig, width = 183, height  = 200, units = "mm")
ggsave(plot = combo, ext_data_sdm_fig_jpg, width = 183, height  = 200, units = "mm", dpi = 300)


# Extended data- robustness of slim simulations ---------------------------
additive_consK_2locus_late <- read_csv(additive_consK_2locus_late_file) %>% 
  mutate(dominance = "additive",
         K_change = "constant",
         loci = 2)

recessive_consK_2locus_late <- read_csv(recessive_consK_2locus_late_file) %>% 
  mutate(dominance = "recessive",
         K_change = "constant",
         loci = 2)

additive_varyK_2locus_late <- read_csv(additive_varyK_2locus_late_file) %>% 
  mutate(dominance = "additive",
         K_change = "varying",
         loci = 2)

recessive_consK_1locus_late <- read_csv(recessive_consK_1locus_late_file) %>% 
  mutate(dominance = "recessive",
         K_change = "constant",
         loci = 1)

selection_key <- tibble(fitness_width = c(0.5446, 0.6562, 0.7805)) %>% 
  mutate(mismatch_penalty = factor(fitness_width, 
                                   levels = c(0.7805, 0.6562, 0.5446),
                                   labels = c("5%", "7%", "10%")))

# Importance of genetic architecture ---------------------------------------
# additive vs recessive vs locus number

for_arch_plot <- rbind(additive_consK_2locus_late, recessive_consK_2locus_late) %>% 
  filter(init_corin == init_ednrb) %>% 
  dplyr::select(-init_ednrb, -freq_ednrb)

arch_plot <- rbind(for_arch_plot, recessive_consK_1locus_late) %>% 
  filter(lambda == 15) %>% 
  rename(p_i = init_corin) %>% 
  mutate(perc_K = N/K) %>% 
  group_by(generation, p_i, fitness_width, dominance, loci) %>% 
  summarise(low95 = quantile(perc_K, c(0.025)),
            up95 = quantile(perc_K, c(0.975)),
            perc_K = mean(perc_K)) %>% 
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(p_i, fitness_width, dominance, loci)) %>% 
  ggplot(aes(x = generation, y = perc_K, ymin = low95, ymax = up95, fill = mismatch_penalty)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2)  +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  geom_ribbon(aes(group = ID), alpha = 0.5) +
  geom_line(aes(color = mismatch_penalty, group = ID)) +
  facet_rep_grid(dominance + loci ~ p_i, labeller = labeller(.rows = label_both, .cols = label_both), scales = "fixed") +
  theme_cowplot() +
  ylab("Percent of carrying capacity, K") +
  ylim(c(0,1.04)) +
  ggtitle("Effect of genetic architecture on evolutionary rescue") +
  theme(axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        plot.title =  element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        strip.text.x  = element_text(size = 5, margin = margin(0.5,0,0.75,0, "mm")),
        strip.text.y = element_text(size = 5, margin = margin(0.75,0,5,0, "mm")),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.title = element_text(size = 5, margin = margin (0,0,0,0, "mm")),
        legend.position = c(0.8, 0.1),
        legend.key.size = unit(3, "mm"))

# Effect of varying carrying capacity -------------------------------------

varyK_plot <- additive_varyK_2locus_late %>% 
  dplyr::select(-min_K, -max_K, -period, -init_K) %>% 
  rbind(mutate(additive_consK_2locus_late, init_dec = NA)) %>% 
  filter(lambda == 15,
         init_corin == init_ednrb) %>% 
  rename(p_i = init_corin) %>% 
  mutate(perc_K = N/K) %>% 
  group_by(generation, p_i, fitness_width, init_dec, K_change) %>% 
  summarise(low95 = quantile(perc_K, c(0.025)),
            up95 = quantile(perc_K, c(0.975)),
            perc_K = mean(perc_K)) %>% 
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(p_i, fitness_width, init_dec, K_change)) %>% 
  ggplot(aes(x = generation, y = perc_K, ymin = low95, ymax = up95, fill = mismatch_penalty)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2)  +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  geom_ribbon(aes(group = ID), alpha = 0.5) +
  geom_line(aes(color = mismatch_penalty, group = ID)) +
  facet_rep_grid(K_change + init_dec ~ p_i, labeller = labeller(.rows = label_both, .cols = label_both), scales = "fixed") +
  theme_cowplot() +
  ylab("Percent of carrying capacity, K") +
  ylim(c(0,1.04)) +
  ggtitle("Effect of population cycling on evolutionary rescue") +
  theme(axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        plot.title =  element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        strip.text.x  = element_text(size = 5, margin = margin(0.5,0,0.75,0, "mm")),
        strip.text.y = element_text(size = 5, margin = margin(0.75,0,5,0, "mm")),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.title = element_text(size = 5, margin = margin (0,0,0,0, "mm")),
        legend.position = c(0.8, 0.1),
        legend.key.size = unit(3, "mm"))

# Combo plot of these two effects -----------------------------------------
combo <- cowplot::plot_grid(arch_plot, varyK_plot,
                            nrow = 2, ncol = 1,
                            labels = c("a", "b"), 
                            label_fontfamily = "sans", 
                            label_fontface = "bold", 
                            label_size = 8,
                            align = "none")


ggsave(plot = combo, filename = arch_varyK_fig_pdf, width = 183, height  = 208, units = "mm")
ggsave(plot = combo, filename = arch_varyK_fig_jpeg, width = 183, height  = 208, units = "mm", dpi = 300)



# Effect of changes in other parameters -----------------------------------

robust <- additive_consK_2locus_late %>% 
  mutate(perc_K = N/K) %>% 
  group_by(generation, init_corin, init_ednrb, fitness_width, lambda) %>% 
  summarise(low95 = quantile(perc_K, c(0.025)),
            up95 = quantile(perc_K, c(0.975)),
            perc_K = mean(perc_K)) %>% 
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(init_corin, init_ednrb, fitness_width, lambda)) %>% 
  ggplot(aes(x = generation, y = perc_K, ymin = low95, ymax = up95, fill = mismatch_penalty)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2)  +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  geom_ribbon(aes(group = ID), alpha = 0.5) +
  geom_line(aes(color = mismatch_penalty, group = ID)) +
  facet_rep_grid(init_ednrb + lambda ~ init_corin, labeller = labeller(.rows = label_both, .cols = label_both), scales = "fixed") +
  theme_cowplot() +
  ylab("Percent of carrying capacity, K") +
  ylim(c(0,1.04)) +
  ggtitle("Effect of initial allele frequencies and reproductive output on evolutionary rescue") +
  theme(axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        plot.title =  element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        strip.text.x  = element_text(size = 5, margin = margin(0.5,0,0.75,0, "mm")),
        strip.text.y = element_text(size = 5, margin = margin(0.75,0,5,0, "mm")),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.title = element_text(size = 5, margin = margin (0,0,0,0, "mm")),
        legend.position = c(0.8, 0.94),
        legend.key.size = unit(3, "mm"))

ggsave(plot = robust, filename = robust_fig_pdf, width = 183, height  = 225, units = "mm")
ggsave(plot = robust, filename = robust_fig_jpeg, width = 183, height  = 225, units = "mm", dpi = 300)

