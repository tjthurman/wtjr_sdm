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
library(rangeBuilder)
library(maptools)
library(geosphere)
library(mapproj)

# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 predicted phenos in rangemap raster

map_file <- args[1]


# Load data -----------------------------------------------------------
#map_file <- "results/pheno/current_predicted_probWhite_SDMrange.tif"
pred.new <- raster(map_file)


state_prov <- ne_states(c("united states of america", "canada"))
countries <- ne_countries(continent = "north america", scale = 10)

# Plot US map -----------------------------------------------------------
pal <- colorRampPalette(colors = c("#69431A", "floralwhite")) # Colors from Mafalda's figs

pdf(file = "results/figures/current_pheno_map.pdf", width = 9.5, height = 7)
plot(
  state_prov,
  xlim = c(-131,-87.5),
  ylim = c(32, 56),
  axes = TRUE,
  col = rgb(133,141,147, maxColorValue = 255),
  border = NA
)

# Plot predicted SDM
plot(pred.new, 
     add = TRUE, 
     legend = F, 
     col = pal(100))
addRasterLegend(pred.new, ramp = pal(100), side = 2,  
                location = c(-127.5, -125.75, 37, 43.5), 
                ncolors = 100, minmax = c(0,1), nTicks = 1)


# Add the state lines back on, but faintly
plot(state_prov, add = T, col = NA, border = "grey61")
# Add countries and coastlines in black
plot(countries, add = T, col = NA, border = "grey10")
plot(coastline10, add = T, col = "grey10", border = "grey10")
box()
maps::map.scale(x = -131, y = 33, ratio = F, cex = .8, relwidth = .11)
dev.off()


# Colorado Insert ---------------------------------------------------------
# Convert raster to df for plotting, subset just to Colorado
pred.new.df <- as(pred.new, "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>% 
  dplyr::filter(x < -102 & x > -109) %>% 
  dplyr::filter(y < 41 & y > 37) %>% 
  dplyr::filter(!is.na(.[1]))
names(pred.new.df) <- c("probWhite", "Long", "Lat")


# Get locations of genomics samples
spec.colo <- read.delim("raw_data/DMNS_spectrometry_PCs.txt", sep = "\t", stringsAsFactors = F) %>% 
  separate(Sample, into = c("museum", "collection", "ID"), remove = F) %>% 
  mutate(colorgroups = cut(PC1, 100)) %>% 
  mutate(colorindex = as.numeric(colorgroups)) %>% 
  mutate(color = pal(100)[colorindex]) 

# Narrow down to unique GPS points, tallying # of samples per GPS
unique_gps <- spec.colo %>% 
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
ggplot() +
  aes(x = Long, y = Lat) +
  coord_quickmap(
    xlim = c(-109, -102),
    ylim = c(37, 41),
    clip = "on", 
    expand = F) +
  geom_tile(data = pred.new.df, aes(fill = probWhite, color = probWhite)) +
  theme(panel.background = element_rect(fill = rgb(133,141,147, maxColorValue = 255), colour = rgb(133,141,147, maxColorValue = 255)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.92, 0.5),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(color = "white", size = 12, face = "bold"),
        legend.title = element_text(color = "white", size = 14, face = "bold")) +
  scale_fill_gradient(low = "#69431A", high  = "floralwhite", breaks = seq(0, 1, length.out = 100), guide = F) +
  scale_color_gradient(low = "#69431A", high  = "floralwhite", breaks = seq(0, 1, length.out = 100), guide = F) +
  geom_point(data  = cluster_plots, aes(size = samples), color = "black", fill = "#F4742A", shape = 21, stroke = 1) +
  scale_size_continuous(range= c(3, 10)) +
  ggsave(filename = "results/figures/colorado.pdf", width = 7, height = 6)


  




# Old colorado plotting ---------------------------------------------------
# Can delete this stuff eventually

# background colors tried
# tomato1 - no good
# slategray1- pretty good, maybe too blue?
# slategray2- pretty good, candidate. But maybe a touch too blue?
# slategray3- pretty good, candidate. But maybe a touch too blue?
# slategray4- pretty good, candidate. Like that a lot
# slategray- pretty good, candidate. Like this a lot too. Think it is the top candidate
# rgb(133,141,147, maxColorValue = 255) A lightened, slightly desautrated slategray. THe new top candidate
# steelblue2 - too blue
# "#82AAD6"- this light blue from Mafalda's figures. I think a little too blue.
# #2F4C9E- the dark blue from Mafalda's figures. I think too blue. 
#F4742A - The orange from Mafalda's figs
# plum2- OK, but i expect it will clash with other figs. 
# plum1- OK, but i expect it will clash with other figs. 
# lightcoral- OK, but I think it would clash. 
# bisque2- Too tan/neutral, I think
# bisque1- Too tan/neutral, I think





pdf(file = "results/figures/colorado.pdf", width = 7, height = 6)
# Start off plot with colorado only
# location only
plot(colorado,
  xlim = c(-109,-102),
  ylim = c(37, 41),
  axes = F,
  col = rgb(133,141,147, maxColorValue = 255),
  border = NA
)
plot(
  state_prov,
  col = rgb(133,141,147, maxColorValue = 255),
  add = T
)
plot(
  colorado,
  col = rgb(133,141,147, maxColorValue = 255),
  add = T
)
plot(pred.new, # or pheno.colorado for colorado only
     add = T,
     col = pal(100),
     legend = F
)
plot(state_prov, add = T, col = NA, border = "grey10")

# Add specc'ed points
# points(x = spec.colo$Long,
#        y = spec.colo$Lat, 
#        pch = 21, col = "black" , bg = spec.colo$color, cex = 3) # colored by PC1
# points(x = spec.colo$Long,
#        y = spec.colo$Lat, 
#        pch = 21, col = "black" , bg = "#F4742A", cex = 2) # or by a different color

points(x = cluster_plots$Long,
       y = cluster_plots$Lat,
       pch = 21, col = "black" , bg = "#F4742A", cex = cluster_plots$samples)
symbols(x = cluster_plots$Long,
       y = cluster_plots$Lat,
       circles = cluster_plots$samples/2, fg = "black" , bg = "#F4742A", add = T)
text(x = cluster_plots$Long,
     y = cluster_plots$Lat,
     labels = cluster_plots$samples)
maps::map.scale(x = -104.5, y = 37.5, ratio = F, col = "black", cex = 1)
box()
dev.off()



# # Get a raster of just colorado
# colorado <- state_prov[state_prov$name == "Colorado",]
# temp <- brick("processed_data/bioclim_30arcsec_for_WTJR_SDM.tif")
# colorado.raster <- rasterize(colorado, temp, field = 1)
# 
# 
# 
# 
# # Mask phenotype predictions down to colorado
# pheno.colorado <- mask(pred.new, colorado.raster)