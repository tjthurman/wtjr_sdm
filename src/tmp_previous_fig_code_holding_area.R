# Horizontal conservation panel --------------------------------------------------
load(conservation_overlap)

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