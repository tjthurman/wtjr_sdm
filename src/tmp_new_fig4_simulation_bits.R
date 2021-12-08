library(tidyverse) 
library(cowplot)
library(lemon)
library(RColorBrewer)
library(ggtext)

selection_key <- tibble(fitness_width = c(0.5446, 0.6562, 0.7805)) %>% 
  mutate(mismatch_penalty = factor(fitness_width, 
                                   levels = c(0.7805, 0.6562, 0.5446),
                                   labels = c("5%", "7%", "10%")))

additive <- read_csv("results/slim_summaries/additive_constantK.csv") %>% 
  mutate(dominance = "additive")

recessive <- read_csv("results/slim_summaries/recessive_constantK.csv") %>% 
  mutate(dominance = "recessive")


# Population trajectories:

bind_rows(additive, recessive) %>% 
  filter(lambda == 15) %>%
  filter(fitness_width == 0.6562) %>% 
  filter(init_corin == init_ednrb) %>% 
  group_by(generation, init_corin, init_ednrb, fitness_width, dominance) %>% 
  summarise(low95 = quantile(N, c(0.025)),
            up95 = quantile(N, c(0.975)),
            N = mean(N)) %>% 
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(init_corin, init_ednrb, dominance)) %>% 
  ggplot(aes(x = generation, y = N, ymin = low95, ymax = up95, fill = as.factor(init_corin))) +
  guides(fill=guide_legend(title="p_i"), color = guide_legend(title = "p_i")) +
  geom_hline(aes(yintercept = 2750), linetype = "dotted") +
  geom_ribbon(aes(group = ID), alpha = 0.5) +
  geom_line(aes(color = as.factor(init_corin), group = ID)) +
  theme_cowplot() +
  theme(legend.position = c(0.1, 0.15)) +
  facet_rep_grid(dominance ~ ., scales = "free")

ggsave("results/figures/tmp/mock_sims_N.pdf", width = 3.4*2, height = 3.35*2, units = "in")

# Heatmap
bind_rows(additive, recessive) %>% 
  filter(lambda == 15) %>%
  filter(init_corin == init_ednrb) %>% 
  filter(generation == 60) %>% 
  group_by(init_corin, init_ednrb, fitness_width, dominance) %>% 
  summarize(deltaN = 2750 - mean(N)) %>% 
  ungroup() %>% 
  left_join(selection_key) %>% 
  ggplot(aes(x = as.factor(init_corin), y = mismatch_penalty, fill = deltaN)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred") +
  facet_grid(dominance ~ .) +
  theme_cowplot()

ggsave("results/figures/tmp/mock_deltaN_heatmap.pdf", width = 3.4*2, height = 3.35*2, units = "in")



# Pop trajectories for Jeff presentation ----------------------------------

bind_rows(additive, recessive) %>% 
  filter(lambda == 15) %>%
  filter(fitness_width == 0.6562) %>% 
  filter(init_corin == init_ednrb) %>% 
  group_by(generation, init_corin, init_ednrb, fitness_width, dominance) %>% 
  summarise(low95 = quantile(N, c(0.025)),
            up95 = quantile(N, c(0.975)),
            N = mean(N)) %>%
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(init_corin, init_ednrb, dominance)) %>% 
  ggplot(aes(x = generation, y = N, ymin = low95, ymax = up95, fill = as.factor(init_corin))) +
  scale_color_viridis_d(begin = 0, end = 0.8) +
  scale_fill_viridis_d(begin = 0, end = 0.8) +
  guides(fill=guide_legend(title="initial<br><i>p</i><sub> brown</sub>"), color = guide_legend(title = "initial<br><i>p</i><sub> brown</sub>")) +
  geom_segment(aes(x = 0, xend = 60, y = 2750, yend = 2750), linetype = "dotted") +
  geom_ribbon(aes(group = ID), alpha = 0.3) +
  geom_line(aes(color = as.factor(init_corin), group = ID), size = 1.5) +
  theme_cowplot() +
  theme(legend.position = c(0.1, 0.15)) +
  facet_rep_grid(dominance ~ .) +
  xlim(c(0,61)) +
  ylab("population size") +
  theme(axis.title = element_text(size = 18),
        legend.title = element_markdown(size = 16),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 16),
        strip.background = element_rect(color = NA, fill = "transparent"))

  

ggsave("results/figures/tmp/sim_pop_traj_for_jeff.pdf", width = 6, height = 7.5, units = "in")



# Plots for figure 4 ------------------------------------------------------

# Panel B using both additive and recessive
bind_rows(additive, recessive) %>% 
  filter(lambda == 15) %>%
  filter(fitness_width == 0.6562) %>% 
  filter(init_corin == init_ednrb) %>% 
  group_by(generation, init_corin, init_ednrb, fitness_width, dominance) %>% 
  summarise(low95 = quantile(N, c(0.025)),
            up95 = quantile(N, c(0.975)),
            N = mean(N)) %>%
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(init_corin, init_ednrb, dominance)) %>% 
  ggplot(aes(x = generation, y = N, ymin = low95, ymax = up95, fill = as.factor(init_corin))) +
  scale_color_viridis_d(begin = 0, end = 0.8) +
  scale_fill_viridis_d(begin = 0, end = 0.8) +
  guides(fill=guide_legend(title="<i>p</i><sub> brown</sub>"), color = guide_legend(title = "<i>p</i><sub> brown</sub>")) +
  geom_segment(aes(x = 0, xend = 60, y = 2750, yend = 2750), linetype = "dotted", size = 0.25) +
  geom_ribbon(aes(group = ID), alpha = 0.3) +
  geom_line(aes(color = as.factor(init_corin), group = ID), size = 0.25) +
  theme_cowplot() +
  theme(legend.position = c(0.07, 0.18)) +
  facet_rep_grid(dominance ~ .) +
  xlim(c(0,61)) +
  #ylab("population size") +
  #ylab("") +
  #xlab("") +
  theme(#axis.title = element_text(size = 5, margin= margin(0,0,0,0, unit = "pt"), lineheight = 0),
        axis.text.y = element_text(size = 4, angle = 90, hjust = 0.5, margin= margin(0,0,0,0, unit = "pt")),
        axis.text.x = element_text(size = 4, margin= margin(0,0,0,0, unit = "pt"), lineheight = 0),
        axis.title = element_blank(),
        axis.ticks.length = unit(0.25, "mm"),
        axis.ticks = element_line(size = 0.25), 
        axis.line = element_line(size = 0.25),
        legend.title = element_markdown(size = 5, margin= margin(0,0,0,0, unit = "pt"),  padding = unit(c(0,0,0,0), "pt")),
        legend.text = element_text(size = 5, margin= margin(0,0,0,0, unit = "pt")),
        legend.key.size = unit(1, "mm"),
        legend.spacing = unit(0, "mm"),
        legend.box.margin = margin(0,0,0,0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.margin = margin(0,0,0,0, unit = "pt"),
        #strip.text = element_text(size = 5, margin= margin(0,0,0,0, unit = "pt")),
        strip.text = element_blank(),
        panel.spacing = unit(0, "mm"),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        plot.margin = margin(0,0,0,0, unit = "pt"),
        plot.background = element_blank(), 
        strip.background = element_blank(),
        strip.placement = "inside")

ggsave("results/figures/sim_pop_trajectories.pdf", width = 31.5, height = 30.5, units = "mm")



# Panel B using only additive
additive %>% 
  filter(lambda == 15) %>%
  filter(fitness_width == 0.6562) %>% 
  filter(init_corin == init_ednrb) %>% 
  group_by(generation, init_corin, init_ednrb, fitness_width) %>% 
  summarise(low95 = quantile(N, c(0.025)),
            up95 = quantile(N, c(0.975)),
            N = mean(N)) %>%
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(init_corin, init_ednrb)) %>% 
  ggplot(aes(x = generation, y = N, ymin = low95, ymax = up95, fill = as.factor(init_corin))) +
  scale_color_viridis_d(begin = 0, end = 0.8) +
  scale_fill_viridis_d(begin = 0, end = 0.8) +
  guides(fill=guide_legend(title="<i>p</i><sub> brown</sub>"), color = guide_legend(title = "<i>p</i><sub> brown</sub>")) +
  geom_segment(aes(x = 0, xend = 60, y = 2750, yend = 2750), linetype = "dotted", size = 0.25) +
  geom_ribbon(aes(group = ID), alpha = 0.3) +
  geom_line(aes(color = as.factor(init_corin), group = ID), size = 0.25) +
  theme_cowplot() +
  theme(legend.position = c(0.07, 0.25)) +
  xlim(c(0,61)) +
  #ylab("population size") +
  #ylab("") +
  #xlab("") +
  theme(#axis.title = element_text(size = 5, margin= margin(0,0,0,0, unit = "pt"), lineheight = 0),
    axis.text.y = element_text(size = 4, angle = 90, hjust = 0.5, margin= margin(0,0,0,0, unit = "pt")),
    axis.text.x = element_text(size = 4, margin= margin(1,0,0,0, unit = "pt"), lineheight = 0),
    axis.title = element_blank(),
    axis.ticks.length = unit(0.25, "mm"),
    axis.ticks = element_line(size = 0.25), 
    axis.line = element_line(size = 0.25),
    legend.title = element_markdown(size = 5, margin= margin(0,0,0,0, unit = "pt"),  padding = unit(c(0,0,0,0), "pt")),
    legend.text = element_text(size = 5, margin= margin(0,0,0,0, unit = "pt")),
    legend.key.size = unit(1.5, "mm"),
    legend.spacing = unit(0, "mm"),
    legend.box.margin = margin(0,0,0,0, unit = "pt"),
    legend.box.spacing = unit(0, "mm"),
    legend.margin = margin(0,0,0,0, unit = "pt"),
    #strip.text = element_text(size = 5, margin= margin(0,0,0,0, unit = "pt")),
    strip.text = element_blank(),
    panel.spacing = unit(0, "mm"),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    plot.margin = margin(0,0,0,0, unit = "pt"),
    plot.background = element_blank(), 
    strip.background = element_blank(),
    strip.placement = "inside")

ggsave("results/figures/sim_pop_trajectories_additive_only.pdf", width = 31.5, height = 30.5, units = "mm")



bind_rows(additive, recessive) %>% 
  filter(lambda == 15) %>%
  filter(fitness_width == 0.6562) %>% 
  filter(init_corin == init_ednrb) %>% 
  group_by(init_corin, init_ednrb, fitness_width, dominance, replicate) %>% 
  mutate(decline = 1 - N/K) %>% 
  slice_max(decline, n = 1) %>% 
  ungroup() %>% 
  group_by(init_corin, init_ednrb, fitness_width, dominance) %>% 
  summarize(mean_max_decline = mean(decline)) %>% 
  ungroup() %>% 
  select(-fitness_width) %>% 
  pivot_wider(names_from = "dominance", values_from = "mean_max_decline") %>% 
  select(-init_ednrb) %>% 
  rename(initial_p_brown = init_corin)

# calculate the percentages lost
bind_rows(additive, recessive) %>% 
  filter(lambda == 15) %>%
  filter(fitness_width == 0.6562) %>% 
  filter(init_corin == init_ednrb) %>% 
  filter(generation == 60) %>% 
  group_by(init_corin, init_ednrb, fitness_width, dominance, replicate) %>% 
  mutate(decline = 1 - N/K) %>% 
  ungroup() %>% 
  group_by(init_corin, init_ednrb, fitness_width, dominance) %>% 
  summarize(mean_max_decline = mean(decline)) %>% 
  ungroup() %>% 
  select(-fitness_width) %>% 
  pivot_wider(names_from = "dominance", values_from = "mean_max_decline") %>% 
  select(-init_ednrb) %>% 
  rename(initial_p_brown = init_corin)



