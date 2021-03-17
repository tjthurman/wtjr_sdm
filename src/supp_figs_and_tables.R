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



# # for running in a pipline
args = commandArgs(trailingOnly=TRUE)
# 1 snow cover glm rdata
# 2 srt glm rdata
# 3 conservation pheno overlap
# 4 current SRT predicted phenos
# 5 future SRT predicted phenos
# 6 current snow cover predicted phenos

# 7 snow cover table out
# 8 snow cover metrics out
# 9 srt table out
# 10 srt metrics out
# 11 broad chisq out
# 12 pheno comparison map out
# 13 glm method compare map out
# 14 percent brown change map out
# 15 discrete current pheno map out

# Inputs
snow_cover_glm <- args[1]
srt_glm <- args[2]
conservation_overlap <- args[3]
current_file <- args[4]
future_file <- args[5]
current_cover_file <- args[6]

# Outputs
snow_cover_table <- args[7]
snow_cover_metrics <- args[8]
srt_table <- args[9]
srt_metrics  <- args[10]
broad_chisq_out <- args[11]
pheno_compare_maps <- args[12]
glm_method_compare_map <- args[13]
percent_brown_change <- args[14]
discrete_current_pheno_map <- args[15]

# For running as a script
# Inputs
# snow_cover_glm <- "results/pheno/current_pheno_glm.RData"
# srt_glm <- "results/pheno/current_pheno_glm_SRT.RData"
# conservation_overlap <- "results/conservation/cons_by_current_color.RData"
# current_file <- "results/pheno/current_predicted_probWhite_SDMrange_SRT.tif"
# future_file <- "results/pheno/future_predicted_probWhite_SDMrange.tif"
# current_cover_file <- "results/pheno/current_predicted_probWhite_SDMrange.tif"
# # outputs
# snow_cover_table <- "results/pheno/glm_table_current_snow_cover.csv"
# snow_cover_metrics <- "results/pheno/glm_metrics_current_snow_cover.csv"
# srt_table <- "results/pheno/glm_table_current_srt.csv"
# srt_metrics <- "results/pheno/glm_metrics_current_srt.csv"
# broad_chisq_out <- "results/conservation/broad_chisq_res.csv"
# pheno_compare_maps <- "results/figures/supplemental/pheno_compare_maps.pdf"
# glm_method_compare_map <- "results/figures/supplemental/model_difference_map.pdf"
# percent_brown_change <- "results/figures/supplemental/percent_brown_change.pdf"
# discrete_current_pheno_map <- "results/figures/supplemental/discrete_current_pheno_map.pdf"

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
future_pheno <- raster(future_file)
current_pheno <- raster(current_file)

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
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255)) +
  geom_tile(data = current_probBrown_df, aes(fill = probBrown, color = probBrown, x = Long, y = Lat)) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, ) +
  geom_sf(data = countries, color = "grey10", fill = NA) +
  geom_sf(data= coast, color = "grey10", fill = NA) +
  coord_sf(
    xlim = c(-125, -89.5),
    ylim = c(32, 49.5),
    clip = "on", 
    expand = F) +
  scale_fill_gradientn(colors = pal(100), guide = F) +
  scale_color_gradientn(colors = pal(100), name = NULL) +
  guides(color = guide_colorbar(label.position = "right", ticks = T,, ticks.colour = "black", frame.colour = "black", frame.linewidth = 1.3, barwidth = 1.5, barheight = 6.75)) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(34, 48, by = 4)) + 
  scale_x_continuous(breaks = seq(-120, -90, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.88, 0.24),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scalebar(x.min = -125, x.max = -89.5,
           y.min = 32, y.max = 49.5, 
           dist  = 200, dist_unit = "km", model = "WGS84", transform = T, anchor = c(x = -92.15, y = 33.3)) +
  ggtitle("Current Prob(Brown)")

map_pheno_future <- ggplot() +
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255)) +
  geom_tile(data = future_probBrown_df, aes(fill = probBrown, color = probBrown, x = Long, y = Lat)) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, ) +
  geom_sf(data = countries, color = "grey10", fill = NA) +
  geom_sf(data= coast, color = "grey10", fill = NA) +
  coord_sf(
    xlim = c(-125, -89.5),
    ylim = c(32, 49.5),
    clip = "on", 
    expand = F) +
  scale_fill_gradientn(colors = pal(100), guide = F) +
  scale_color_gradientn(colors = pal(100), name = NULL) +
  guides(color = guide_colorbar(label.position = "right", ticks = T,, ticks.colour = "black", frame.colour = "black", frame.linewidth = 1.3, barwidth = 1.5, barheight = 6.75)) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(34, 48, by = 4)) + 
  scale_x_continuous(breaks = seq(-120, -90, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.88, 0.24),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scalebar(x.min = -125, x.max = -89.5,
           y.min = 32, y.max = 49.5, 
           dist  = 200, dist_unit = "km", model = "WGS84", transform = T, anchor = c(x = -92.15, y = 33.3)) +
  ggtitle("Future Prob(Brown)")

combo <- grid.arrange(map_pheno_current,
             map_pheno_future,
             nrow = 2) 
ggsave(plot = combo, pheno_compare_maps, width = 9.7, height  = 14, units = "in")


# Difference in SRT vs. snow cover models ---------------------------------
current_snow_cover <- raster(current_cover_file)
current_snow_cover_masked <- mask(crop(current_snow_cover, current_pheno), current_pheno)
  
# Get model difference, do it in terms of probBrown
model_difference <- (1- current_snow_cover_masked) - (1 - current_pheno)

model_difference_current_probBrown <- model_difference %>% 
  as(., "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1])) 
names(model_difference_current_probBrown) <- c("difference", "Long", "Lat")

diff_pal <- colorRampPalette(colors = c("#762a83", "#f7f7f7", "#1b7837"))

map_model_difference <- ggplot() +
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255)) +
  geom_tile(data = model_difference_current_probBrown, aes(fill = difference, color = difference, x = Long, y = Lat)) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, ) +
  geom_sf(data = countries, color = "grey10", fill = NA) +
  geom_sf(data= coast, color = "grey10", fill = NA) +
  coord_sf(
    xlim = c(-125, -89.5),
    ylim = c(32, 49.5),
    clip = "on", 
    expand = F) +
  scale_fill_gradientn(colors = diff_pal(100), guide = F) +
  scale_color_gradientn(colors = diff_pal(100), name = NULL) +
  guides(color = guide_colorbar(label.position = "right", ticks = T,, ticks.colour = "black", frame.colour = "black", frame.linewidth = 1.3, barwidth = 1.5, barheight = 6.75)) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(34, 48, by = 4)) + 
  scale_x_continuous(breaks = seq(-120, -90, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.88, 0.24),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scalebar(x.min = -125, x.max = -89.5,
           y.min = 32, y.max = 49.5, 
           dist  = 200, dist_unit = "km", model = "WGS84", transform = T, anchor = c(x = -92.15, y = 33.3)) +
  ggtitle("Difference in current Prob(Brown) between snow cover and SRT models")

ggsave(plot = map_model_difference, glm_method_compare_map, width = 9.7, height  = 7, units = "in")


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

area_compare %>% 
  pivot_longer(current:future, names_to = "time_period", values_to = "perc.brown") %>% 
  ggplot(aes(x = prob.brown.thresh, y = perc.brown, color = time_period)) +
  geom_line() +
  xlab("Threshold of Prob(brown)") +
  ylab("Percent of WTJR range above threshold (\"brown\")") +
  scale_color_discrete(type = c(alpha("purple", 1), alpha("darkgreen", 1))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  ggsave(percent_brown_change, height = 4, width = 4, units = "in")


# Map of current, discrete phenos -----------------------------------------
pred.new <- raster(current_cover_file)

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
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(35, 55, by = 5)) + 
  scale_x_continuous(breaks = seq(-130, -70, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.08, 0.36),
        legend.key = element_rect(colour = "black"),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")) +
  ggsn::scalebar(x.min = -132, x.max = -80.5,
                 y.min = 32, y.max = 56.5, 
                 dist  = 300, dist_unit = "km", model = "WGS84", transform = T, anchor = c(x = -122.48, y = 35))


ggsave(discrete_pheno, filename = discrete_current_pheno_map, width = 9.7, height = 7)