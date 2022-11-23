
# Load packages -----------------------------------------------------------
library(rnaturalearth)
library(rnaturalearthhires)
library(raster)
library(ggsn)
library(geosphere)
library(readxl)
library(ggokabeito)
library(tidyverse)


# Set up map stuff --------------------------------------------------------
state_prov <- ne_states(c("united states of america", "canada"), returnclass = "sf")
countries <- ne_countries(continent = "north america", scale = 10, returnclass = "sf")
coast <- ne_coastline(scale = 10, returnclass = "sf")



# Read in IUCN range ------------------------------------------------------
iucn_range <- sf::read_sf("raw_data/redlist_species_data_79e8f518-a14c-4d0e-9640-5bea84d7c1b8/data_0.shp")

# Get SDM occurrences -----------------------------------------------------
in_csv <- "processed_data/thin/wtjr_occ_0km_thin1.csv"
sdm_occurrences <- read.csv(in_csv, stringsAsFactors = F) %>%
  dplyr::select(longitude = roundlon, latitude = roundlat) %>% 
  mutate(`Sample type` = "SDM- presence")



# Get SDM rangemap --------------------------------------------------------
range_raster <- raster("results/sdm/sdm_rangemap_best_sens95.grd") %>% 
  aggregate(., fact = 8)


z <- rasterToPolygons(range_raster, dissolve = T)

# Do some modification of the polygones to be able to plot with ggplot
# in order to plot polygons, first fortify the data
z@data$id <- rownames(z@data)
# create a data.frame from our spatial object
zdata <- fortify(z, region = "id")
# merge the "fortified" data with the data from our spatial object
range_df <- merge(zdata, z@data,
            by = "id") %>% 
  mutate(name = ifelse(!hole, "SDM range", "hole")) %>% 
  arrange(hole, order) %>% 
  mutate(order = 1:nrow(.))

holes <- range_df %>% 
  filter(hole)



# get locations for all genetic samples -----------------------------------
sample_info <- read_csv("raw_data/genetics/samples_info.csv")

# Data on which samples went into the 74 for GWAS
disambig <- read_xlsx("../wtjr_genomics/data/WTJR_74lowcovsamples_code_disamb.xlsx") 


genetic_samples <- sample_info %>% 
  select(sample, latitude, longitude) %>% 
  mutate(`Sample type` = ifelse(sample %in% str_replace(disambig$SampleID, "_", "-"), "WGS- GWAS", "WGS- other"))



# Combine it all
to_plot <- bind_rows(sdm_occurrences,
                     genetic_samples)


# Plot --------------------------------------------------------------------
# Pick better colors

state_fill <- "grey91"
sdm_border_col <- "grey31"
sdm_fill <- "grey80"
iucn_color <- palette_okabe_ito()[9]  


ggplot() +
  geom_sf(data = state_prov, color = "grey51", fill = state_fill, size = 0.25) +
  geom_polygon(data = range_df, aes(x = long, y = lat, group = group, fill = name), color = sdm_border_col) + 
  geom_polygon(data = holes, aes(x = long, y = lat, group = group), fill = state_fill, color = sdm_border_col) + 
  scale_fill_manual(values = c("SDM range" = sdm_fill,
                               "hole" = state_fill),
                    limits = "SDM range") +
  geom_point(data = to_plot, aes(x = longitude, y = latitude, 
                                 color = `Sample type`, 
                                 size = `Sample type`)) +
  geom_sf(data = state_prov, color = "grey51", fill = NA, size = 0.25) +
  geom_sf(data = countries, color = "grey10", fill = NA, size = 0.25) +
  geom_sf(data = coast, color = "grey10", fill = NA, size = 0.25) +
  geom_sf(data = iucn_range, aes(lty = "IUCN range"), color = iucn_color, fill = NA, size = 1) +
  coord_sf(
    xlim = c(-132, -87.8),
    ylim = c(32, 60.5),
    clip = "on", 
    expand = F)  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.84, 0.84),
        legend.box.background  = element_rect(fill = "grey98", color = "black", size = 0.8),
        legend.background = element_blank(),
        legend.spacing = unit(-0.5, "mm"),
        legend.box.margin = margin(1,1,1,1, unit = "pt"),
        legend.margin = margin(1,3,1,3, unit = "pt"),
        legend.key = element_blank(),
        axis.title = element_blank(), 
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_continuous(breaks = seq(-130, -89.2, by = 10)) +
  scale_size_manual(values = c("SDM- presence" = 0.75,
                               "WGS- GWAS" = 2,
                               "WGS- other" = 2)) +
  scale_color_manual(values = c("SDM- presence" = palette_okabe_ito()[1],
                                "WGS- GWAS" = palette_okabe_ito()[5],
                                "WGS- other" = palette_okabe_ito()[3])) +
  guides(fill = guide_legend(title = NULL, title.position = "bottom", override.aes = list(color = sdm_border_col)),
         color = guide_legend(title.theme = element_text(size = 10)),
         linetype = guide_legend(title = NULL, title.position = "bottom")) +
  ggsn::scalebar(x.min = -132, x.max = -90,
                 y.min = 32, y.max = 60.5, 
                 dist  = 300, dist_unit = "km", model = "WGS84", transform = T, 
                 anchor = c(x = -123.5, y = 33.5), st.size = 10/.pt, border.size = 0.25)
ggsave("results/figures/supplemental/sampling_map.pdf", width = 6.25, height = 6.25, units = "in")
ggsave("results/figures/supplemental/sampling_map.jpeg", width = 6.25, height = 6.25, units = "in", dpi = 320)


