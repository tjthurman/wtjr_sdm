


# OLD
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




# Old main map pre-insert -------------------------------------------------

us <- ggplot() +
  geom_sf(data = state_prov, color = "grey61", fill = rgb(133,141,147, maxColorValue = 255)) +
  geom_tile(data = change.df, aes(fill = change_probBrown, color = change_probBrown, x = Long, y = Lat)) +
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
  guides(color = guide_colorbar(label.position = "right", ticks = T, ticks.colour = "black", frame.colour = "black", frame.linewidth = 1.3, barwidth = 1.5, barheight = 6.75)) +
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

