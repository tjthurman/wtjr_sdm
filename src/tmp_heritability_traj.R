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

# Load results from the additive simulations, done with constant and new VE
additive_consK_2locus_early <- read_csv("results/slim_summaries_80gens/additive_constantK_2locus_early.csv") %>% 
  mutate(dominance = "additive",
         K_change = "constant",
         loci = 2,
         Ve = "recalculated",
         scenario = "additive_constantK_2locus")



additive_consK_2locus_early_consVE <- read_csv("results/slim_summaries_80gens_consVE/additive_constantK_2locus_early.csv") %>% 
  mutate(dominance = "additive",
         K_change = "constant",
         loci = 2,
         Ve = "constant",
         scenario = "additive_constantK_2locus")

recessive_consK_2locus_early_consVE <- read_csv("results/slim_summaries_80gens_consVE/recessive_constantK_2locus_early.csv") %>% 
  mutate(dominance = "recessive",
         K_change = "constant",
         loci = 2,
         Ve = "constant",
         scenario = "recessive_constantK_2locus")

recessive_consK_1locus_early_consVE <- read_csv("results/slim_summaries_80gens_consVE/recessive_constantK_1locus_early.csv") %>% 
  mutate(dominance = "recessive",
         K_change = "constant",
         loci = 1,
         Ve = "constant",
         scenario = "recessive_constantK_1locus")

additive_varyK_2locus_early_consVE <- read_csv("results/slim_summaries_80gens_consVE/additive_varyK_2locus_early.csv") %>% 
  mutate(dominance = "additive",
         K_change = "varying",
         loci = 2,
         Ve = "constant",
         scenario = "additive_varyK_2locus")


# Plot heritability trajectories ------------------------------------------
# First, compare how heritability changes across the old way and the new way

old_vs_new <- rbind(additive_consK_2locus_early, additive_consK_2locus_early_consVE) %>% 
  filter(lambda == 15) %>% 
  filter(fitness_width == 0.6562) %>% 
  group_by(generation, init_corin, init_ednrb, Ve) %>% 
  summarise(low95 = quantile(obs_H2, c(0.025)),
            up95 = quantile(obs_H2, c(0.975)),
            obs_H2 = mean(obs_H2)) %>% 
  ungroup() %>%
  mutate(ID = paste(init_corin, init_ednrb, Ve)) %>% 
  mutate(Ve = as.factor(Ve))

old_vs_new_plot <- old_vs_new %>% 
  ggplot(aes(x = generation, y = obs_H2, ymin = low95, ymax = up95, fill = Ve)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2)  +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  geom_hline(aes(yintercept = 0.64), linetype = "dashed") +
  geom_ribbon(aes(group = ID), alpha = 0.5) +
  geom_line(aes(color = Ve, group = ID)) +
  facet_rep_grid(init_corin ~ init_ednrb , labeller = labeller(.rows = label_both, .cols = label_both), scales = "free") +
  theme_cowplot() +
  ylab("Heritability (Va/Vp)") +
  ggtitle("Effect of Ve modelling on heritability") +
  theme(axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        plot.title =  element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        strip.text.x  = element_text(size = 5, margin = margin(0.5,0,0.75,0, "mm")),
        strip.text.y = element_text(size = 5, margin = margin(0.75,0,5,0, "mm")),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.title = element_text(size = 5, margin = margin (0,0,0,0, "mm")),
        legend.key.size = unit(3, "mm"))

ggsave(old_vs_new_plot, filename = "results/figures/tmp/additive_2locus_heritability_old_vs_new.pdf", width = 30, height = 30, units = "cm")


# Then, in the new way, compare how the heritability changes with changes in allele frequency
h2_vs_allele_freqs <- additive_consK_2locus_early_consVE %>% 
  filter(lambda == 15) %>% 
  filter(fitness_width == 0.6562) %>%
  pivot_longer(cols = c(freq_corin, freq_ednrb, obs_H2)) %>% 
  group_by(generation, init_corin, init_ednrb, Ve, name) %>% 
  summarise(low95 = quantile(value, c(0.025)),
            up95 = quantile(value, c(0.975)),
            value = mean(value)) %>% 
  ungroup() %>% 
  mutate(ID = paste(init_corin, init_ednrb, Ve, name)) %>% 
  mutate(Ve = as.factor(Ve))
  

h2_vs_allele_freqs_plot <- h2_vs_allele_freqs %>% 
  ggplot(aes(x = generation, y = value, ymin = low95, ymax = up95, fill = name)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2)  +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  geom_hline(aes(yintercept = 0.64), linetype = "dashed") +
  geom_ribbon(aes(group = ID), alpha = 0.5) +
  geom_line(aes(color = name, group = ID)) +
  facet_rep_grid(init_corin ~ init_ednrb , labeller = labeller(.rows = label_both, .cols = label_both), scales = "free") +
  theme_cowplot() +
  ylab("Heritability (Va/Vp)") +
  ggtitle("Trajectories of allele frequencies and heritability") +
  theme(axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        plot.title =  element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        strip.text.x  = element_text(size = 5, margin = margin(0.5,0,0.75,0, "mm")),
        strip.text.y = element_text(size = 5, margin = margin(0.75,0,5,0, "mm")),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.title = element_text(size = 5, margin = margin (0,0,0,0, "mm")),
        legend.key.size = unit(3, "mm"))

ggsave(h2_vs_allele_freqs_plot, filename = "results/figures/tmp/additive_2locus_heritability_vs_allele_freq.pdf", width = 30, height = 30, units = "cm")



# Finally, look across the 4 genetic architectures ------------------------

a <- additive_varyK_2locus_early_consVE %>% 
  dplyr::select(-min_K, -max_K, -period, -init_K) %>% 
  dplyr::select(init_dec, init_ednrb, freq_ednrb, everything())

b <- additive_consK_2locus_early_consVE %>% 
  mutate(init_dec = NA) %>% 
  dplyr::select(init_dec, init_ednrb, freq_ednrb, everything())

c <- recessive_consK_1locus_early_consVE %>% 
  mutate(init_dec = NA, freq_ednrb = NA, init_ednrb = NA) %>% 
  dplyr::select(init_dec, init_ednrb, freq_ednrb, everything())

d <- recessive_consK_2locus_early_consVE %>% 
  mutate(init_dec = NA) %>% 
  dplyr::select(init_dec, init_ednrb, freq_ednrb, everything())


h2_across_scenarios <- rbind(a,b,c,d) %>% 
  filter(lambda == 15) %>% 
  filter(fitness_width == 0.6562) %>%
  group_by(generation, init_corin, init_ednrb, init_dec, scenario) %>% 
  summarise(low95 = quantile(obs_H2, c(0.025), na.rm = T),
            up95 = quantile(obs_H2, c(0.975), na.rm = T),
            value = mean(obs_H2, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(ID = paste(init_corin, init_ednrb, init_dec, scenario)) 

h2_across_scenarios_plot <- h2_across_scenarios %>% 
  ggplot(aes(x = generation, y = value, ymin = low95, ymax = up95, fill = scenario)) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2)  +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  geom_hline(aes(yintercept = 0.64), linetype = "dashed") +
  geom_ribbon(aes(group = ID), alpha = 0.5) +
  geom_line(aes(color = scenario, group = ID)) +
  facet_rep_grid(init_corin ~ init_ednrb , labeller = labeller(.rows = label_both, .cols = label_both), scales = "free") +
  theme_cowplot() +
  ylab("Heritability (Va/Vp)") +
  ggtitle("Effect of genetic architecture on heritability") +
  theme(axis.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        axis.title = element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        plot.title =  element_text(size = 7, margin = margin(0,0,0,0, "mm")),
        strip.text.x  = element_text(size = 5, margin = margin(0.5,0,0.75,0, "mm")),
        strip.text.y = element_text(size = 5, margin = margin(0.75,0,5,0, "mm")),
        legend.text = element_text(size = 5, margin = margin(0,0,0,0, "mm")),
        legend.title = element_text(size = 5, margin = margin (0,0,0,0, "mm")),
        legend.key.size = unit(3, "mm"))
ggsave(h2_across_scenarios_plot, filename = "results/figures/tmp/h2_across_scenarios.pdf", width = 30, height = 30, units = "cm")




