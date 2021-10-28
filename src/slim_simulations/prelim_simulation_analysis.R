library(tidyverse)
library(cowplot)

late <- NULL
for (file in list.files("slim_results/additive_constantK/", pattern = "late.csv", full.names = T)) {
  # parse filename for sim parameters
  elements <- str_split(basename(file), pattern = "_")[[1]]
  sim_gens <- as.numeric(str_remove(elements[1], "gens"))
  init_Corin <- as.numeric(str_remove(elements[4], "iCorin"))
  init_EDNRB <- as.numeric(str_remove(elements[5], "iEDNRB"))
  init_opt <- as.numeric(str_remove(elements[6], "initOpt"))
  final_opt <- as.numeric(str_remove(elements[7], "finalOpt"))
  fitness_width <- as.numeric(str_remove(elements[8], "sel"))
  lambda  <- as.numeric(str_remove(elements[10], "lambda"))
  replicate <- as.numeric(str_remove(elements[11], "rep"))
  
  # Import CSV and add parameters
  one_sim_res <- read.csv(file) %>% 
    mutate(sim_gens = sim_gens,
           init_corin = init_Corin,
           init_ednrb = init_EDNRB,
           fitness_width = fitness_width,
           lambda = lambda,
           replicate = replicate)
  late <- rbind(late, one_sim_res)
}




# Population size through time
late %>% 
  mutate(ID = paste(init_corin, init_ednrb, fitness_width, lambda, replicate)) %>% 
  filter(init_corin == init_ednrb) %>% 
  ggplot(aes(x = generation, y = N, color = as.factor(fitness_width))) +
  geom_line(aes(group = ID), alpha = 0.4) +
  geom_smooth() +
  ylim(c(0, 5250)) +
  theme_cowplot() +
  geom_hline(aes(yintercept = 5000)) +
  facet_grid(init_corin ~ lambda)


opt_pheno_df <- late %>% 
  filter(init_corin == init_ednrb,
         init_corin == 0.01,
         fitness_width == 0.6562,
         lambda == 4,
         replicate ==1) %>% 
  select(generation, opt_pheno)

# Average phenotype through time, compared to optimum
late %>% 
  filter(init_corin == init_ednrb,
         init_corin == 0.01,
         fitness_width == 0.6562,
         lambda == 4) %>% 
  ggplot(aes(x = generation, y = mean_pheno)) +
  geom_line(aes(group = replicate), alpha = 0.5) +
  geom_smooth() +
  geom_line(aes(x = generation, y = opt_pheno), linetype = "dashed", data = opt_pheno_df)  +
  theme_cowplot() 




# Max population decrease




# Some gut check metrics:
# percent survival, age structure, etc

