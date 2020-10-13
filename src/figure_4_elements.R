####################
# Project: WTJR SDMs
# Author: Timothy Thurman
# Purpose: Create figure elements for figure 4
# Date Created: Tue Apr 21 15:05:30 2020
####################

# Load packages -----------------------------------------------------------
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(raster)
library(rangeBuilder)
library(maptools)
library(geosphere)
library(mapproj)
library(sf)
library(ggsn)
library(cowplot)
library(lemon)

# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 future predicted phenos
# 2 conservation/pheno overlap

map_file <- args[1]
conservation_file <- args[2]


future_file <- "results/pheno/future_predicted_probWhite_SDMrange.tif"
current_file <- "results/pheno/current_predicted_probWhite_SDMrange_SRT.tif"

# Load data ---------------------------------------------------------------
future_pheno <- raster(future_file)
current_pheno <- raster(current_file)
state_prov <- ne_states(c("united states of america", "canada"), returnclass = "sf")
countries <- ne_countries(continent = "north america", scale = 10, returnclass = "sf")
coast <- ne_coastline(scale = 10, returnclass = "sf")


state_prov_base <- ne_states(c("united states of america", "canada"))
countries_base <- ne_countries(continent = "north america", scale = 10)


load("results/conservation/cons_by_current_color.RData")

# Plot US map -----------------------------------------------------------
pal <- colorRampPalette(colors = c("white", "#F4742A")) # orange
pal <- colorRampPalette(colors = c("white", "#664A9E")) # purple
pal <- colorRampPalette(colors = c("#323333", "#C061FF")) # mafalda suggested, tested
pal <- colorRampPalette(colors = c("#190064", "#FFF700")) # mafalda suggested, tested.
pal <- colorRampPalette(colors = c("#190064", "#FFBB00")) # mafalda suggested, tested
pal <- colorRampPalette(colors = c("#33164D", "#00A9E0")) # mafalda suggested, tested.
pal <- colorRampPalette(colors = c(rgb(41, 57, 113, maxColorValue = 255), # mills 2018
                                   rgb(248, 236, 63, maxColorValue = 255),
                                   rgb(176, 85, 59, maxColorValue = 255)))
  

difference <- current_pheno - future_pheno


# A ggplot way
future.df <- as(future_pheno, "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1]))
names(future.df) <- c("probWhite", "Long", "Lat")

change.df <- as(difference, "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1]))
names(change.df) <- c("change_probWhite", "Long", "Lat")

us <- ggplot() +
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255)) +
  geom_tile(data = change.df, aes(fill = change_probWhite, color = change_probWhite, x = Long, y = Lat)) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, ) +
  geom_sf(data = countries, color = "grey10", fill = NA) +
  geom_sf(data= coast, color = "grey10", fill = NA) +
  coord_sf(
    xlim = c(-125, -89.5),
    ylim = c(32, 49.5),
    clip = "on", 
    expand = F) +
  scale_fill_gradientn(colors = pal(100), guide = F) +
  scale_color_gradientn(colors = pal(100), name = NULL, breaks = c(0, 0.4, 0.8)) +
  guides(color = guide_colorbar(label.position = "right", ticks = F, frame.colour = "black", frame.linewidth = 1.3, barwidth = 1.5, barheight = 6.75)) +
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
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")) +
  scalebar(x.min = -125, x.max = -89.5,
           y.min = 32, y.max = 49.5, 
           dist  = 200, dist_unit = "km", model = "WGS84", transform = T, anchor = c(x = -92.15, y = 33.3)) 

us 

us +
  ggsave(filename = "results/figures/future_pheno_map.pdf", width = 9.7, height = 7)

us +
  ggsave(filename = "results/figures/future_pheno_map.png", width = 9.7, height = 7, dpi = 72)


# Conservation status bar plot --------------------------------------------
broad$conservation <- factor(broad$conservation, 
       levels = c("extirpated", "broad.extirp", "local.extirp", "poss.decline", "pres.stable"),
       labels = c("extirpated", "broad\nextirpations", "local\nextirpations",
                  "possible\ndeclines", "presumed\nstable"))
broad %>% 
  pivot_longer(brown:mixed, names_to = "pheno", values_to = "area") %>% 
  ggplot(aes(x = conservation, y = area/1000, fill = pheno)) +
  geom_col(color = "grey30") +
  theme_cowplot() +
  scale_fill_manual(values = c(
    "brown" = "#69431A",
    "mixed" = "#B39D83", # midpoint of our palette
    "white" = "floralwhite"
  ), name = "phenotype") +
  xlab("conservation status") +
  ylab(expression(paste("area (k", m^2, "), thousands"))) +
  theme(legend.position = c(0.035, 0.85),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 20),
        axis.text.y.left  = element_text(size = 18, angle = 90, hjust = 0.5),
        legend.key.size = unit(1.5, "cm")) +
  ggsave(filename = "results/figures/cons_barplot_test.pdf", width = 9, height = 8, unit = "in")

broad %>% 
  pivot_longer(brown:mixed, names_to = "pheno", values_to = "area") %>% 
  ggplot(aes(x = pheno, y = area/1000, fill = conservation)) +
  geom_col(color = "grey30") +
  scale_fill_manual(
    values = c(
      "extirpated" = "#EC4741",
      "broad\nextirpations" = "#F8B732",
      "local\nextirpations" = "#F2F031",
      "possible\ndeclines" = "#ABD344",
      "presumed\nstable" = "#44854B"
    ), name = "conservation\nstatus"
  ) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  theme_cowplot() +
  xlab("phenotype") +
  ylab(expression(paste("area (k", m^2, "), thousands"))) +
  theme(legend.position = c(0.035, 0.75),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 20),
        axis.text.y.left  = element_text(size = 18, angle = 90, hjust = 0.5),
        legend.key.size = unit(1.5, "cm")) +
  ggsave(filename = "results/figures/pheno_barplot_test.pdf", width = 8, height = 8, unit = "in")



# A third, artier option --------------------------------------------------


broad
  

