
# Load libraries ----------------------------------------------------------
library(tidyverse) 
library(cowplot)
library(lemon)
library(RColorBrewer)
library(ggtext)
library(ggrepel)
library(broom)


# Load selection key ------------------------------------------------------
selection_key <- tibble(fitness_width = c(0.5446, 0.6562, 0.7805)) %>% 
  mutate(mismatch_penalty = factor(fitness_width, 
                                   levels = c(0.7805, 0.6562, 0.5446),
                                   labels = c("5%", "7%", "10%")))

# Load results from the 4 sets of simulations -----------------------------
additive_consK_2locus_late <- read_csv("results/slim_summaries_80gens_consVE/additive_constantK_2locus_late.csv") %>% 
  mutate(dominance = "additive",
         K_change = "constant",
         loci = 2)

recessive_consK_2locus_late <- read_csv("results/slim_summaries_80gens_consVE/recessive_constantK_2locus_late.csv") %>% 
  mutate(dominance = "recessive",
         K_change = "constant",
         loci = 2)

additive_varyK_2locus_late <- read_csv("results/slim_summaries_80gens_consVE/additive_varyK_2locus_late.csv") %>% 
  mutate(dominance = "additive",
         K_change = "varying",
         loci = 2)

recessive_consK_1locus_late <- read_csv("results/slim_summaries_80gens_consVE/recessive_constantK_1locus_late.csv") %>% 
  mutate(dominance = "recessive",
         K_change = "constant",
         loci = 1)

# Important of genetic architecture ---------------------------------------
# additive vs recessive vs locus number

for_arch_plot <- rbind(additive_consK_2locus_late, recessive_consK_2locus_late) %>% 
  filter(init_corin == init_ednrb) %>% 
  dplyr::select(-init_ednrb, -freq_ednrb)

arch_plot <- rbind(for_arch_plot, recessive_consK_1locus_late) %>% 
  filter(lambda == 15) %>% 
  rename(p_i = init_corin) %>% 
  mutate(perc_K = N/K) %>% 
  group_by(generation, p_i, fitness_width, dominance, loci) %>% 
  summarise(low95 = quantile(perc_K, c(0.025)),
            up95 = quantile(perc_K, c(0.975)),
            perc_K = mean(perc_K)) %>% 
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(p_i, fitness_width, dominance, loci)) %>% 
  ggplot(aes(x = generation, y = perc_K, ymin = low95, ymax = up95, fill = mismatch_penalty)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2)  +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  geom_ribbon(aes(group = ID), alpha = 0.5) +
  geom_line(aes(color = mismatch_penalty, group = ID)) +
  facet_rep_grid(dominance + loci ~ p_i, labeller = labeller(.rows = label_both, .cols = label_both), scales = "fixed") +
  theme_cowplot() +
  ylab("Percent of carrying capacity, K") +
  ylim(c(0,1.04)) +
  ggtitle("Effect of genetic architecture on evolutionary rescue") +
  theme(axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        plot.title =  element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        strip.text.x  = element_text(size = 5, margin = margin(0.5,0,0.75,0, "mm")),
        strip.text.y = element_text(size = 5, margin = margin(0.75,0,5,0, "mm")),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.title = element_text(size = 5, margin = margin (0,0,0,0, "mm")),
        legend.position = c(0.8, 0.1),
        legend.key.size = unit(3, "mm"))

# Effect of varying carrying capacity -------------------------------------

varyK_plot <- additive_varyK_2locus_late %>% 
  dplyr::select(-min_K, -max_K, -period, -init_K) %>% 
  rbind(mutate(additive_consK_2locus_late, init_dec = NA)) %>% 
  filter(lambda == 15,
         init_corin == init_ednrb) %>% 
  rename(p_i = init_corin) %>% 
  mutate(perc_K = N/K) %>% 
  group_by(generation, p_i, fitness_width, init_dec, K_change) %>% 
  summarise(low95 = quantile(perc_K, c(0.025)),
            up95 = quantile(perc_K, c(0.975)),
            perc_K = mean(perc_K)) %>% 
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(p_i, fitness_width, init_dec, K_change)) %>% 
  ggplot(aes(x = generation, y = perc_K, ymin = low95, ymax = up95, fill = mismatch_penalty)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2)  +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  geom_ribbon(aes(group = ID), alpha = 0.5) +
  geom_line(aes(color = mismatch_penalty, group = ID)) +
  facet_rep_grid(K_change + init_dec ~ p_i, labeller = labeller(.rows = label_both, .cols = label_both), scales = "fixed") +
  theme_cowplot() +
  ylab("Percent of carrying capacity, K") +
  ylim(c(0,1.04)) +
  ggtitle("Effect of population cycling on evolutionary rescue") +
  theme(axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        plot.title =  element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        strip.text.x  = element_text(size = 5, margin = margin(0.5,0,0.75,0, "mm")),
        strip.text.y = element_text(size = 5, margin = margin(0.75,0,5,0, "mm")),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.title = element_text(size = 5, margin = margin (0,0,0,0, "mm")),
        legend.position = c(0.8, 0.1),
        legend.key.size = unit(3, "mm"))

# Combo plot of these two effects -----------------------------------------
combo <- cowplot::plot_grid(arch_plot, varyK_plot,
                            nrow = 2, ncol = 1,
                            labels = c("a", "b"), 
                            label_fontfamily = "sans", 
                            label_fontface = "bold", 
                            label_size = 8,
                            align = "none")


ggsave(plot = combo, filename = "results/figures/supplemental/extended_data_arch_varyK_sim_res_80gens_consVE.pdf", width = 183, height  = 208, units = "mm")
ggsave(plot = combo, filename = "results/figures/supplemental/extended_data_arch_varyK_sim_res_80gens_consVE.jpeg", width = 183, height  = 208, units = "mm", dpi = 300)



# Effect of changes in other parameters -----------------------------------

robust <- additive_consK_2locus_late %>% 
  mutate(perc_K = N/K) %>% 
  group_by(generation, init_corin, init_ednrb, fitness_width, lambda) %>% 
  summarise(low95 = quantile(perc_K, c(0.025)),
            up95 = quantile(perc_K, c(0.975)),
            perc_K = mean(perc_K)) %>% 
  ungroup() %>% 
  left_join(selection_key) %>%
  mutate(ID = paste(init_corin, init_ednrb, fitness_width, lambda)) %>% 
  ggplot(aes(x = generation, y = perc_K, ymin = low95, ymax = up95, fill = mismatch_penalty)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2)  +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  geom_ribbon(aes(group = ID), alpha = 0.5) +
  geom_line(aes(color = mismatch_penalty, group = ID)) +
  facet_rep_grid(init_ednrb + lambda ~ init_corin, labeller = labeller(.rows = label_both, .cols = label_both), scales = "fixed") +
  theme_cowplot() +
  ylab("Percent of carrying capacity, K") +
  ylim(c(0,1.04)) +
  ggtitle("Effect of initial allele frequencies and reproductive output on evolutionary rescue") +
  theme(axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        plot.title =  element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        strip.text.x  = element_text(size = 5, margin = margin(0.5,0,0.75,0, "mm")),
        strip.text.y = element_text(size = 5, margin = margin(0.75,0,5,0, "mm")),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.title = element_text(size = 5, margin = margin (0,0,0,0, "mm")),
        legend.position = c(0.8, 0.94),
        legend.key.size = unit(3, "mm"))

ggsave(plot = robust, filename = "results/figures/supplemental/extended_data_sim_robust_80gens_consVE.pdf", width = 183, height  = 225, units = "mm")

ggsave(plot = robust, filename = "results/figures/supplemental/extended_data_sim_robust_80gens_consVE.jpeg", width = 183, height  = 225, units = "mm", dpi = 300)




