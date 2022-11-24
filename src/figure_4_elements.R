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
library(lemon)
library(ggtext)


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
slim_res_additive <- args[3]
slim_res_recessive <- args[4]
map_out_pdf <- args[5]
map_out_png <- args[6]
insert_plot_out <- args[7]
pop_traj_plot_out <- args[8]
perc_brown_csv_out <- args[9]

# for running as a script
# current_file <- "results/pheno/current_predicted_probWhite_SDMrange_SRT.tif"
# future_file <- "results/pheno/future_predicted_probWhite_SDMrange.tif"
# slim_res_additive <- "results/slim_summaries/additive_constantK_2locus_late.csv"
# slim_res_recessive <- "results/slim_summaries/recessive_constantK_2locus_late.csv"
# map_out_pdf <- "results/figures/pheno_change_map_55mm.pdf"
# map_out_png <- "results/figures/pheno_change_map_55mm.png"
# insert_plot_out <- "results/figures/density_probBrown_insert_55mm.pdf"
# pop_traj_plot_out <- "results/figures/sim_pop_trajectories_55mm.pdf"
# perc_brown_csv_out <- "results/pheno/percent_brown_by_time.csv"

# Load data ---------------------------------------------------------------
future_pheno <- raster(future_file)
current_pheno <- raster(current_file)
state_prov <- ne_states(c("united states of america", "canada"), returnclass = "sf")
countries <- ne_countries(continent = "north america", scale = 10, returnclass = "sf")
coast <- ne_coastline(scale = 10, returnclass = "sf")

# Plot US map -----------------------------------------------------------
pal <- colorRampPalette(colors = c(rgb(41, 57, 113, maxColorValue = 255), # mills 2018  USING THIS AT THE MOMENT
                                   rgb(248, 236, 63,  alpha = 0.8, maxColorValue = 255),
                                   rgb(212, 59, 15, maxColorValue = 255)))

# Full resolution
difference <- (1- future_pheno) - (1 - current_pheno) 

# Lower res for plotting
difference_lowres <- aggregate(difference, fact = 4)

change.df <- as(difference_lowres, "SpatialPixelsDataFrame") %>% 
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
  scale_x_continuous(breaks = seq(-125, -105, by = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.095, 0.34),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_blank(),
        axis.ticks = element_line(size = 0.35, color = "black"),
        axis.ticks.length = unit(0.5, "mm"),
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


ggsave(us, filename = map_out_png, width = 54.5, height = 34, units = "mm", dpi = 300)
ggsave(us, filename = map_out_pdf, width = 54.5, height = 34, units = "mm")

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
        plot.margin = margin(0.08, 0.1, 0, 0.04, "cm"))

ggsave(insert, filename = insert_plot_out, height = 13.8, width = 16.5, units = "mm")



# Simulation trajectories -------------------------------------------------
additive <- read_csv(slim_res_additive) %>% 
  mutate(dominance = "additive")
recessive <- read_csv(slim_res_recessive) %>% 
  mutate(dominance = "recessive")

selection_key <- tibble(fitness_width = c(0.5446, 0.6562, 0.7805)) %>% 
  mutate(mismatch_penalty = factor(fitness_width, 
                                   levels = c(0.7805, 0.6562, 0.5446),
                                   labels = c("5%", "7%", "10%")))


# Panel B using both additive and recessive
trajectories <- bind_rows(additive, recessive) %>% 
  filter(lambda == 15) %>%
  filter(fitness_width == 0.6562) %>% 
  filter(init_corin == init_ednrb) %>% 
  mutate(perc_K = N/K) %>% 
  group_by(generation, init_corin, init_ednrb, fitness_width, dominance) %>% 
  summarise(low95 = quantile(perc_K, c(0.025)),
            up95 = quantile(perc_K, c(0.975)),
            perc_K = mean(perc_K)) %>% 
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(init_corin, init_ednrb, dominance)) %>% 
  ggplot(aes(x = generation, y = perc_K, ymin = low95, ymax = up95, fill = as.factor(init_corin))) +
  scale_color_viridis_d(begin = 0, end = 0.8) +
  scale_fill_viridis_d(begin = 0, end = 0.8) +
  guides(fill=guide_legend(title="<i>p</i><sub> i</sub>"), color = guide_legend(title = "<i>p</i><sub> i</sub>")) +
  geom_segment(aes(x = 0, xend = 81, y = 1, yend = 1), linetype = "dotted", size = 0.25) +
  geom_ribbon(aes(group = ID), alpha = 0.3) +
  geom_line(aes(color = as.factor(init_corin), group = ID), size = 0.25) +
  theme_cowplot() +
  theme(legend.position = c(0.03, 0.35)) +
  facet_rep_grid(. ~ dominance, scales = "fixed", repeat.tick.labels = T) +
  xlim(c(0,81)) +
  ylim(c(0.4, 1.04)) +
  theme(
    axis.text.y = element_text(size = 5, angle = 90, hjust = 0.5, margin= margin(0,0,0,0, unit = "pt")),
    axis.text.x = element_text(size = 5, margin= margin(1,0,0,0, unit = "pt"), lineheight = 0),
    axis.title = element_blank(),
    axis.ticks.length = unit(0.35, "mm"),
    axis.ticks = element_line(size = 0.35), 
    axis.line = element_line(size = 0.35),
    legend.title = element_markdown(size = 5, margin= margin(0,0,0,0, unit = "pt"),  padding = unit(c(0,0,0,0), "pt")),
    legend.text = element_text(size = 5, margin= margin(0,0,0,0, unit = "pt")),
    legend.key.size = unit(1.5, "mm"),
    legend.spacing = unit(0, "mm"),
    legend.box.margin = margin(0,0,0,0, unit = "pt"),
    legend.box.spacing = unit(0, "mm"),
    legend.margin = margin(0,0,0,0, unit = "pt"),
    strip.text = element_blank(),
    panel.spacing = unit(0, "mm"),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    plot.margin = margin(0,0,0,0, unit = "pt"),
    plot.background = element_blank(), 
    strip.background = element_blank(),
    strip.placement = "inside")

ggsave(trajectories, filename = pop_traj_plot_out, width = 52.5, height = 19, units = "mm")


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
