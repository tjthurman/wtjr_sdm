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
library(ggsn)
library(cowplot)


# Get arguments -----------------------------------------------------------
# # for running in a pipline
args = commandArgs(trailingOnly=TRUE)
# 1 current SRT predicted phenos
# 2 future SRT predicted phenos
# 3 conservation overlap Rdata
# 4 change map pdf
# 5 change map png
# 6 conservation plot pdf
# 7 insert plot pdf
# 8 percent brown out csv

current_file <- args[1]
future_file <- args[2]
conservation_overlap <- args[3]
map_out_pdf <- args[4]
map_out_png <- args[5]
cons_plot_out <- args[6]
insert_plot_out <- args[7]
perc_brown_csv_out <- args[8]

# for running as a script
# current_file <- "results/pheno/current_predicted_probWhite_SDMrange_SRT.tif"
# future_file <- "results/pheno/future_predicted_probWhite_SDMrange.tif"
# conservation_overlap <- "results/conservation/cons_by_current_color.RData"
# map_out_pdf <- "results/figures/pheno_change_map_83mm.pdf"
# map_out_png <- "results/figures/pheno_change_map_83mm.png"
# cons_plot_out <- "results/figures/horizontal_consv.pdf"
# insert_plot_out <- "results/figures/density_probBrown_insert.pdf"
# perc_brown_csv_out <- "results/pheno/percent_brown_by_time.csv"

# Load data ---------------------------------------------------------------
future_pheno <- raster(future_file)
current_pheno <- raster(current_file)
state_prov <- ne_states(c("united states of america", "canada"), returnclass = "sf")
countries <- ne_countries(continent = "north america", scale = 10, returnclass = "sf")
coast <- ne_coastline(scale = 10, returnclass = "sf")

load(conservation_overlap)

# Plot US map -----------------------------------------------------------
pal <- colorRampPalette(colors = c(rgb(41, 57, 113, maxColorValue = 255), # mills 2018  USING THIS AT THE MOMENT
                                   rgb(248, 236, 63,  alpha = 0.8, maxColorValue = 255),
                                   rgb(212, 59, 15, maxColorValue = 255)))


difference <- (1- future_pheno) - (1 - current_pheno)

change.df <- as(difference, "SpatialPixelsDataFrame") %>% 
  as.data.frame(.) %>%
  dplyr::filter(!is.na(.[1]))
names(change.df) <- c("change_probBrown", "Long", "Lat")

us <- ggplot() +
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255), size = 0.25) +
  geom_tile(data = change.df, aes(fill = change_probBrown, color = change_probBrown, x = Long, y = Lat)) +
  geom_sf(data = state_prov, color = "grey61", fill = NA, size = 0.25) +
  geom_sf(data = countries, color = "grey10", fill = NA, size = 0.25) +
  geom_sf(data= coast, color = "grey10", fill = NA, size = 0.25) +
  coord_sf(
    xlim = c(-129.75, -89.5),
    ylim = c(32, 49.5),
    clip = "on", 
    expand = F) +
  scale_fill_gradientn(colors = pal(100), guide = F) +
  scale_color_gradientn(colors = pal(100), name = NULL, breaks = c(0, 0.4, 0.8)) +
  guides(color = guide_colorbar(label.position = "left", ticks = T, ticks.colour = "black", frame.colour = "black", 
                                frame.linewidth = 0.75,
                                barwidth = unit(1.75, units = "mm"),
                                barheight = unit(8, units = "mm"))) +
  scale_y_continuous(breaks = seq(34, 48, by = 4)) + 
  scale_x_continuous(breaks = seq(-125, -95, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.095, 0.34),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_blank(),
        axis.ticks = element_line(size = 0.75, color = "black"),
        axis.text = element_text(size = 5, color = "black"),
        axis.title = element_blank(),
        axis.text.y.left = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 5),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", size = 0.75),
        plot.margin = unit(c(0.05, 0.05, 0, 0), "cm"), 
        panel.spacing = unit(0, "cm"),
        plot.background = element_blank()) +
  scalebar(x.min = -129.75, x.max = -89.5,
           y.min = 32, y.max = 49.5, 
           dist  = 350, dist_unit = "km", model = "WGS84", transform = T, anchor = c(x = -121, y = 33.4),
           st.size = 5/.pt, border.size = 0.25, st.dist = 0.0325) 


ggsave(us, filename = map_out_png, width = 57, height = 34.34, units = "mm", dpi = 72)


ggsave(us, filename = map_out_pdf, width = 57, height = 34.34, units = "mm")






# Horizontal conservation panel --------------------------------------------------
broad$conservation <- factor(broad$conservation, 
                             levels = rev(c("extirpated", "broad.extirp", "local.extirp", "poss.decline", "pres.stable")),
                             labels = rev(c("extirpated", "broad\nextirpations", "local\nextirpations",
                                        "possible\ndeclines", "presumed\nstable")))

color_pal <- c("black", "grey30", "grey60", "grey80", "grey95") # greyscale

broad %>% 
  pivot_longer(brown:mixed, names_to = "pheno", values_to = "area") %>% 
  mutate(pheno = factor(.$pheno, levels = c("white", "mixed", "brown"))) %>% 
  ggplot(aes(x = pheno, y = area/1000, fill = conservation)) +
  geom_col(color = "black", position = position_fill(), width = 0.5, size = 0.15) +
  scale_fill_manual(
    values = c(
      "extirpated" = color_pal[1],
      "broad\nextirpations" = color_pal[2],
      "local\nextirpations" = color_pal[3],
      "possible\ndeclines" = color_pal[4],
      "presumed\nstable" = color_pal[5]
    ), name = "conservation\nstatus", guide = F
  ) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  coord_flip() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  ggsave(cons_plot_out, width = 4.25, height = 2, units = "cm")



# Insert --------------------------------------------------------
# Density of Prob white across years:
cur_pheno_vec <- values(current_pheno) %>% 
  .[!is.na(.)]
future_pheno_vec <- values(future_pheno) %>% 
  .[!is.na(.)]

probBrown_df <- data.frame(time = c(rep("current", times = length(cur_pheno_vec)),
                                    rep("future", times = length(future_pheno_vec))),
                           probWhite = c(cur_pheno_vec, future_pheno_vec)) %>% 
  mutate(probBrown = 1 - probWhite)

insert <- probBrown_df %>% 
  ggplot(aes(x = probBrown, color = time, fill = time)) +
  geom_density(aes(y = ..scaled..), size = 0.25) +
  scale_fill_discrete(type = c(alpha("grey60", 0.5), alpha("#2F4C9E", 0.5))) +
  scale_color_discrete(type = c(alpha("grey30", 1), alpha("darkblue", 0.75))) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0), name = "P(brown)", breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0,1.01), expand = c(0, 0), breaks = c(0,1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.5, 0.8),
        legend.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.y = element_text(color = "black", size = 5, margin = margin(0,0,0,0, "mm")),
        axis.text.x = element_text(color = "black", size = 5, margin = margin(0,0,0,0, "mm"), vjust = 2.5),
        axis.ticks = element_blank(),
        axis.line = element_line(color = "black", size = 0.25),
        axis.title.x = element_text(color = "black", size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 5, color = "black", margin = margin(0,0,0,0, "mm")),
        legend.key.size = unit(5/.pt, "mm"),
        legend.key = element_blank(),
        legend.spacing = unit(0, "mm"),
        legend.margin = margin(0,0,0,0, "mm"), 
        # panel.spacing = unit(0.05, "mm"),
        #panel.border = element_blank(),
        plot.margin = margin(0.08, 0.1, 0, 0.04, "cm"))


ggsave(insert, filename = insert_plot_out, height = 13.8, width = 17.5, units = "mm")



# Stats on proportion brown --------------------------------------------------------
current_for_brown_area <- current_pheno
values(current_for_brown_area) <- ifelse(values(current_for_brown_area) <= 0.2, 1, NA)

future_for_brown_area <- future_pheno
values(future_for_brown_area) <- ifelse(values(future_for_brown_area) <= 0.2, 1, NA)


area_brown_current <- sum(values(area(current_for_brown_area, na.rm = T)), na.rm = T)
area_brown_future <- sum(values(area(future_for_brown_area, na.rm = T)), na.rm = T)

us_only_range_area <- sum(values(area(current_pheno, na.rm = T)), na.rm = T)

current_perc_brown <- area_brown_current/us_only_range_area
future_perc_brown <- area_brown_future/us_only_range_area

data.frame(perc.brown.current = current_perc_brown,
           perc.brown.future = future_perc_brown) %>%
    mutate(fold.change = perc.brown.future/perc.brown.current) %>%
    write.csv(., file = perc_brown_csv_out, row.names = F)
