
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
# Get samples from the bamlist
seqed_samples <- read_delim(file = "../wtjr_genomics/data/all_wtjr_sample_bamlist.txt", delim = "\t", col_names = "bam_file") %>% 
  mutate(sample = basename(bam_file)) %>% 
  mutate(sample = str_remove(sample, ".realigned.sorted.bam")) %>% 
  mutate(sample = str_remove(sample, "_realigned_sorted.bam"))


# My data
tim <- read_csv("../wtjr_genomics/data/labwork/possible_samples_to_extracts.csv") %>% 
  select(sample = ID_museum, year, latitude, longitude, winter_pheno, state, source = source.x) %>% 
  mutate(museum = year < 2008) %>% 
  mutate(museum = ifelse(is.na(museum), T, museum)) %>% 
  mutate(batch = "tim")

# Fix the LTW_TWO_4140 sample
tim$sample[tim$sample == "Mafalda_LTW_WYO_4140"] <- "LTW_COL_4140"
tim$year[tim$sample == "LTW_COL_4140"] <- NA
tim$latitude[tim$sample == "LTW_COL_4140"] <- 39.04814
tim$longitude[tim$sample == "LTW_COL_4140"] <- -105.56139
tim$winter_pheno[tim$sample == "LTW_COL_4140"] <- NA
tim$state[tim$sample == "LTW_COL_4140"] <- "CO"
tim$source[tim$sample == "LTW_COL_4140"] <- "LTW_COL_4140"


# Mafalda's data
mafalda <- read_excel("../wtjr_genomics/data/mf_thesis_table_S3_3.xlsx") %>% 
  filter(Species == "Lepus townsendii") %>%
  mutate(museum = ifelse(Origin == "Field", F, T),
         year = NA, 
         winter_pheno = NA,
         sample = str_replace(`Sample ID`, "_", "-"),
         batch = "mafalda") %>% 
  separate(`Location (Population)`, into = c("country", "state"), sep = "/") %>% 
  select(sample, year, latitude = `Lat.`, longitude = `Long.`,  winter_pheno, source = Origin, state, museum, batch) %>% 
  bind_rows(., read_csv("../wtjr_genomics/data/mafalda_UCM_samples.csv"))

# Data on which samples went into the 74 for GWAS
disambig <- read_xlsx("../wtjr_genomics/data/WTJR_74lowcovsamples_code_disamb.xlsx") 


genetic_samples <- seqed_samples %>% 
  left_join(bind_rows(tim, mafalda)) %>% 
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
  
ggplot() +
  geom_sf(data = state_prov, color = "grey51", fill = state_fill, size = 0.25) +
  geom_polygon(data = range_df, aes(x = long, y = lat, group = group, fill = name), color = sdm_border_col) + 
  geom_polygon(data = holes, aes(x = long, y = lat, group = group), fill = state_fill, color = sdm_border_col) + 
  geom_point(data = to_plot, aes(x = longitude, y = latitude, color = `Sample type`), size = 1.5) +
  geom_sf(data = state_prov, color = "grey51", fill = NA, size = 0.25) +
  geom_sf(data = countries, color = "grey10", fill = NA, size = 0.25) +
  geom_sf(data = coast, color = "grey10", fill = NA, size = 0.25) +
  coord_sf(
    xlim = c(-132, -89.2),
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
  scale_fill_manual(values = c("SDM range" = sdm_fill,
                               "hole" = state_fill),
                    limits = "SDM range") +
  scale_color_okabe_ito() +
  guides(fill = guide_legend(title = NULL, title.position = "bottom", override.aes = list(color = sdm_border_col)),
         color = guide_legend(title.theme = element_text(size = 10))) +
  ggsn::scalebar(x.min = -132, x.max = -90,
                 y.min = 32, y.max = 60.5, 
                 dist  = 300, dist_unit = "km", model = "WGS84", transform = T, 
                 anchor = c(x = -123.5, y = 33.5), st.size = 10/.pt, border.size = 0.25)
ggsave("../wtjr_sdm/results/prelim_sampling_map.pdf")


