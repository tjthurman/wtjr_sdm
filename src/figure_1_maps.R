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


# For tunning as script
# map_file <- "results/pheno/current_predicted_probWhite_SDMrange.tif"
# colo_samples <- "raw_data/sample_coordinates_74individuals.txt"
# us_pdf <- "results/figures/current_pheno_map.pdf"
# us_png <- "results/figures/current_pheno_map.png"
# colo_pdf <- "results/figures/colorado.pdf"
# colo_png <- "results/figures/colorado.png"



# Load data -----------------------------------------------------------
pred.new <- raster(map_file)


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
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255)) +
  geom_tile(data = pred.pheno.df, aes(fill = probBrown, color = probBrown, x = Long, y = Lat)) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, ) +
  geom_sf(data = countries, color = "grey10", fill = NA) +
  geom_sf(data= coast, color = "grey10", fill = NA) +
  coord_sf(
    xlim = c(-132, -80.5),
    ylim = c(32, 56.5),
    clip = "on", 
    expand = F) +
  scale_fill_gradientn(colors = pal(100), guide = F) +
  scale_color_gradientn(colors = pal(100), name = NULL, breaks = c(0, 0.5, 1)) +
  guides(color = guide_colorbar(label.position = "left", ticks = T, ticks.colour = "black", frame.colour = "black", frame.linewidth = 1.3, barwidth = 1.5, barheight = 7)) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = seq(35, 55, by = 5)) + 
  scale_x_continuous(breaks = seq(-130, -100, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.105, 0.36),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")) +
  ggsn::scalebar(x.min = -132, x.max = -80.5,
           y.min = 32, y.max = 56.5, 
           dist  = 300, dist_unit = "km", model = "WGS84", transform = T, anchor = c(x = -122.48, y = 35))



ggsave(us, filename = us_pdf, width = 9.7, height = 7)
ggsave(us, filename = us_png, width = 9.7, height = 7, dpi = 72)

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

# Then, plot. Will use ggplot so I can scale circles by sample size
colorado <- ggplot() +
  aes(x = Long, y = Lat) +
  coord_quickmap(
    xlim = c(-109, -102),
    ylim = c(37, 41),
    clip = "on", 
    expand = F) +
  geom_tile(data = pred.pheno.df, aes(fill = probBrown, color = probBrown)) +
  theme(panel.background = element_rect(fill = rgb(133,141,147, maxColorValue = 255), colour = rgb(133,141,147, maxColorValue = 255)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.89, 0.27),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(color = "white", size = 21, face = "bold"),
        legend.title = element_text(color = "white", size = 23, face = "bold", vjust = 0, margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
        legend.spacing.y = unit(1.5, "pt")) +
  scale_fill_gradientn(colors = pal(100), guide = F) +
  scale_color_gradientn(colors = pal(100), guide = F) +
  geom_point(data  = cluster_plots, aes(size = samples), color = "black", fill = "#F4742A", shape = 21, stroke = 1) +
  scale_size_continuous(range= c(6, 14)) +
  ggsn::scalebar(x.min = -109, x.max = -102,
           y.min = 37, y.max = 41, 
           dist  = 50, dist_unit = "km", model = "WGS84", transform = T, anchor = c(x = -102.6, y = 37.25), st.color = "white", st.size = 6, st.dist = 0.025) 


ggsave(colorado, filename = colo_pdf, width = 7, height = 6)

ggsave(colorado, filename = colo_png, width = 7, height = 6, dpi = 72)

  


