####################
# Project:
# Author: Timothy Thurman
# Purpose:
# Date Created: Tue Apr 21 15:05:30 2020
####################

# Load packages -----------------------------------------------------------
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(raster)
library(rangeBuilder)
library(maptools)

# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 predicted phenos in rangemap raster

map_file <- args[1]


# Load data -----------------------------------------------------------
map_file <- "results/pheno/current_predicted_probWhite_SDMrange.tif"
pred.new <- raster(map_file)


state_prov <- ne_states(c("united states of america", "canada"))
countries <- ne_countries(continent = "north america", scale = 10)
# PLot US map -----------------------------------------------------------
pal <- colorRampPalette(colors = c("#69431A", "floralwhite")) # Colors from Mafalda's figs

pdf(file = "results/figures/current_pheno_map.pdf", width = 9.5, height = 7)
plot(
  state_prov,
  xlim = c(-131,-87.5),
  ylim = c(32, 56),
  axes = TRUE,
  col = "slategray",
  border = NA
)

# colors tried
# tomato1 - no good
# slategray1- pretty good, maybe too blue?
# slategray2- pretty good, candidate. But maybe a touch too blue?
# slategray3- pretty good, candidate. But maybe a touch too blue?
# slategray4- pretty good, candidate. Like that a lot
# slategray- pretty good, candidate. Like this a lot too. Think it is the top candidate
# steelblue2 - too blue
# "#82AAD6"- this light blue from Mafalda's figures. I think a little too blue.
# #2F4C9E- the dark blue from Mafalda's figures. I think too blue. 
#F4742A - The orange from Mafalda's figs
# plum2- OK, but i expect it will clash with other figs. 
# plum1- OK, but i expect it will clash with other figs. 
# lightcoral- OK, but I think it would clash. 
# bisque2- Too tan/neutral, I think
# bisque1- Too tan/neutral, I think


# Plot predicted SDM
plot(pred.new, 
     add = TRUE, 
     legend = F, 
     col = pal(100))
addRasterLegend(pred.new, ramp = pal(100), side = 2,  
                location = c(-128.5, -126, 35, 47), 
                ncolors = 100, minmax = c(0,1), nTicks = 3)


# Add the state lines back on, but faintly
plot(state_prov, add = T, col = NA, border = "grey61")
# Add countries and coastlines in black
plot(countries, add = T, col = NA, border = "grey10")
plot(coastline10, add = T, col = "grey10", border = "grey10")
# Maybe add a box showing the colorado inset??
# plot(colorado, add = T, col = NA, border = "grey5")
box()
maps::map.scale(x = -132.8, y = 33.8, ratio = F)
dev.off()


# Colorado insert

# Colorado Insert ---------------------------------------------------------
# Get a raster of just colorado
colorado <- state_prov[state_prov$name == "Colorado",]
temp <- brick("processed_data/bioclim_30arcsec_for_WTJR_SDM.tif")
colorado.raster <- rasterize(colorado, temp, field = 1)


# Mask phenotype predictions down to colorado
pheno.colorado <- mask(pred.new, colorado.raster)


# Get spec-ed samples that are within colorado
# colorize their PC scores with the same colors as
# the predicted phenotypes
spec.colo <- read.delim("raw_data/DMNS_spectrometry_PCs.txt", sep = "\t", stringsAsFactors = F) %>% 
  separate(Sample, into = c("museum", "collection", "ID"), remove = F) %>% 
  mutate(colorgroups = cut(PC1, 100)) %>% 
  mutate(colorindex = as.numeric(colorgroups)) %>% 
  mutate(color = pal(100)[colorindex]) %>% 
  filter(as.numeric(Lat) > 37) %>% 
  filter(as.numeric(Lat) < 41) %>% 
  filter(as.numeric(Long) < -102) %>% 
  filter(as.numeric(Long) > -109)

eps(file = "results/figures/colorado.pdf", width = 7, height = 7)
# Start off plot with colorado only
# location only
plot(colorado,
  xlim = c(-109,-102),
  ylim = c(37, 41),
  axes = F,
  col = "slategray",
  border = NA
)
plot(
  state_prov,
  col = "slategray",
  add = T
)
plot(
  colorado,
  col = "slategray",
  add = T
)
plot(pheno.colorado, # or pheno.colorado for colorado only
     add = T,
     col = pal(100),
     legend = F
)
plot(state_prov, add = T, col = NA, border = "grey10")

# Add specc'ed points
points(x = spec.colo$Long,
       y = spec.colo$Lat, 
       pch = 21, col = "black" , bg = spec.colo$color, cex = 3) # colored by PC1
points(x = spec.colo$Long,
       y = spec.colo$Lat, 
       pch = 21, col = "black" , bg = "#664A9E", cex = 3) # or by a different color
maps::map.scale(x = -108.7, y = 37.25, ratio = F, col = "black")
box()
dev.off()