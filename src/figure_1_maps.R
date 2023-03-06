####################
# Project: WTJR SDMs
# Author: Timothy Thurman
# Purpose: Create figure elements for figure 1
# Date Created: Tue Apr 21 15:05:30 2020
####################

# Load packages -----------------------------------------------------------
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(raster)
library(ggsn)
library(geosphere)
library(ggrastr)

# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 predicted phenos in rangemap raster
# 2 coordinates of samples in colorado
# 3 US pdf file
# 4 US png file 
# 5 colorado pdf file
# 6 colorago png file

map_file <- args[1]
colo_samples <- args[2]
us_pdf <- args[3]
us_png <- args[4]
colo_pdf <- args[5]
colo_png <- args[6]


# For running as script
# map_file <- "results/pheno/current_predicted_probWhite_SDMrange.tif"
# colo_samples <- "raw_data/sample_coordinates_74individuals.txt"
# us_pdf <- "results/figures/current_pheno_map_120mm.pdf"
# us_png <- "results/figures/current_pheno_map_120mm.png"
# colo_pdf <- "results/figures/colorado_120mm.pdf"
# colo_png <- "results/figures/colorado_120mm.png"

# Load data -----------------------------------------------------------
pred.new <- raster(map_file) %>% 
  aggregate(., fact = 4)


state_prov <- ne_states(c("united states of america", "canada"), returnclass = "sf")
countries <- ne_countries(continent = "north america", scale = 10, returnclass = "sf")
coast <- ne_coastline(scale = 10, returnclass = "sf")

# Plot US map -----------------------------------------------------------
pal <- colorRampPalette(colors = rev(c("#69431A", "floralwhite"))) # Colors from Mafalda's figs


pred.pheno.df <- as(pred.new, "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1])) %>% 
  mutate(probBrown = 1 - current_predicted_probWhite_SDMrange)
names(pred.pheno.df) <- c("probWhite", "Long", "Lat", "probBrown")

us <- ggplot() +
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255), size = 0.25) +
  ggrastr::rasterise(geom_tile(data = pred.pheno.df, aes(fill = probBrown, color = probBrown, x = Long, y = Lat)), dpi = 400) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, size = 0.25) +
  geom_sf(data = countries, color = "grey10", fill = NA, size = 0.25) +
  geom_sf(data= coast, color = "grey10", fill = NA, size = 0.25) +
  coord_sf(
    xlim = c(-132, -80.5),
    ylim = c(32, 56.5),
    clip = "on", 
    expand = F) +
  scale_fill_gradientn(colors = pal(100), guide = F) +
  scale_color_gradientn(colors = pal(100), name = NULL, breaks = c(0, 0.5, 1)) +
  guides(color = guide_colorbar(label.position = "left", 
                                ticks = T, 
                                ticks.colour = "black", 
                                frame.colour = "black", 
                                frame.linewidth = 1, 
                                barwidth = unit(3.5, "mm"), 
                                barheight = unit(14, "mm"),
                                label.hjust = 1)) +
  scale_y_continuous(breaks = seq(35, 55, by = 5)) + 
  scale_x_continuous(breaks = seq(-130, -100, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.09, 0.36),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.text.align = 1, 
        legend.box.margin = margin(0,0,0,0, unit = "pt"),
        legend.margin = margin(0,0,0,0, unit = "pt"),
        legend.key = element_blank(),
        legend.text = element_text(size = 7, margin = margin(0,0,0,0, unit = "pt"), hjust = 1),
        axis.ticks = element_line(size = 1),
        axis.title = element_blank(), 
        axis.text = element_text(size = 7, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  ggsn::scalebar(x.min = -132, x.max = -80.5,
           y.min = 32, y.max = 56.5, 
           dist  = 300, dist_unit = "km", model = "WGS84", transform = T, anchor = c(x = -123.5, y = 33.5), st.size = 7/.pt, border.size = 0.25)

ggsave(us, filename = us_pdf, width = 115, height = 74, units = "mm")
ggsave(us, filename = us_png, width = 115, height = 74, units = "mm", dpi = 300)


# Colorado Insert ---------------------------------------------------------
# Get locations of genomics samples
colo.samples <- read.delim(colo_samples)

# Narrow down to unique GPS points, tallying # of samples per GPS
unique_gps <- colo.samples %>% 
  mutate(roundLat = round(Lat, digits = 5),
         roundLong = round(Long, digits = 5)) %>% 
  group_by(roundLat, roundLong) %>% 
  tally() %>% 
  arrange(roundLat, roundLong) %>% 
  ungroup()

# For next round of clustering, will get the distance of each location to the previous point
unique_gps$dist_to_next <- NA
for (row in 2:dim(unique_gps)[1]) {
  unique_gps$dist_to_next[row] <- distVincentyEllipsoid(p1 = c(unique_gps$roundLong[row - 1],
                                                               unique_gps$roundLat[row - 1]),
                                                        p2 = c(unique_gps$roundLong[row],
                                                               unique_gps$roundLat[row]))/1000
}

# Then cluster them based on that
# Can choose different limits, will try 15km for now. 
dist_lim <- 15
cluster <- 1
unique_gps$cluster <- NA
for (row in 1:dim(unique_gps)[1]) {
  if (row == 1) {
    unique_gps$cluster[row] <- 1
  }  else if (unique_gps$dist_to_next[row] < dist_lim) {
    unique_gps$cluster[row] <- cluster
  } else {
    cluster <- cluster + 1
    unique_gps$cluster[row] <- cluster
  }
}


cluster_plots <- unique_gps %>% 
  group_by(cluster) %>% 
  summarize(Lat = mean(roundLat),
            Long = mean(roundLong),
            samples = sum(n))


# Need to figure out how to match legend key size with the displayed size
# Then, plot. Will use ggplot so I can scale circles by sample size
colorado <- ggplot() +
  aes(x = Long, y = Lat) +
  coord_quickmap(
    xlim = c(-109, -102),
    ylim = c(37, 41.02),
    clip = "on", 
    expand = F) +
  rasterize(geom_tile(data = pred.pheno.df, aes(fill = probBrown, color = probBrown)), dpi = 400) +
  scale_fill_gradientn(colors = pal(100), guide = F) +
  scale_color_gradientn(colors = pal(100), guide = F) +
  geom_point(data  = cluster_plots, aes(size = samples), color = "black", fill = "#F4742A", shape = 21, stroke = 0.25) +
  scale_size_continuous(range= c(1, 4)) +
  theme(panel.background = element_rect(fill = rgb(133,141,147, maxColorValue = 255), colour = rgb(133,141,147, maxColorValue = 255)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.86, 0.245),
        legend.background = element_blank(),
        legend.margin = margin(0,0,0,0, "pt"),
        legend.key = element_blank(),
        legend.key.size = unit(0.2, "mm"),
        legend.text = element_text(color = "white", size = 7, face = "bold", margin = margin(0,0,0,0, unit = "mm")),
        legend.title = element_text(color = "white", size = 7, face = "bold", vjust = 0, margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
        legend.spacing = unit(.03, "mm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0.05, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  labs(x = NULL, y = NULL)

ggsave(colorado, filename = colo_pdf, width = 43.5, height = 32, unit = "mm")
ggsave(colorado, filename = colo_png, width = 43.5, height = 32, unit = "mm", dpi = 300)

  