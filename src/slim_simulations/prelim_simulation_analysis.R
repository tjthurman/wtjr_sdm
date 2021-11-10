library(tidyverse)




# Selection key for plotting ----------------------------------------------
selection_key <- tibble(fitness_width = c(0.5446, 0.6562, 0.7805)) %>% 
  mutate(mismatch_penalty = factor(fitness_width, 
                                   levels = c(0.7805, 0.6562, 0.5446),
                                   labels = c("5%", "7%", "10%")))

# Scenario 1 --------------------------------------------------------------
# Complete mismatch, no SGV
late_mismatch_noSGV <- NULL
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
  
  if (fitness_width == 0.65) {
    next
  }
  
  if (init_Corin == 0 & init_EDNRB == 0 & init_opt == 1 & final_opt == 1) {
    # Import CSV and add parameters
    one_sim_res <- read.csv(file) %>%
      mutate(sim_gens = sim_gens,
             init_corin = init_Corin,
             init_ednrb = init_EDNRB,
             init_opt = init_opt,
             final_opt = final_opt,
             fitness_width = fitness_width,
             lambda = lambda,
             replicate = replicate)
    late_mismatch_noSGV <- rbind(late_mismatch_noSGV, one_sim_res)
  }
}


write_csv(late_mismatch_noSGV, path = "results/tmp/slim_summaries/complete_mismatch_no_SGV.csv", col_names = T)



# Moving optimum no SGV ---------------------------------------------------
late_moveOpt_noSGV <- NULL
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
  
  if (fitness_width == 0.65) {
    next
  }
  
  if (init_Corin == 0 & init_EDNRB == 0 & init_opt != 1) {
    # Import CSV and add parameters
    one_sim_res <- read.csv(file) %>%
      mutate(sim_gens = sim_gens,
             init_corin = init_Corin,
             init_ednrb = init_EDNRB,
             init_opt = init_opt,
             final_opt = final_opt,
             fitness_width = fitness_width,
             lambda = lambda,
             replicate = replicate)
    late_moveOpt_noSGV <- rbind(late_moveOpt_noSGV, one_sim_res)
  }
}


write_csv(late_moveOpt_noSGV, path = "results/tmp/slim_summaries/moving_opt_no_SGV.csv", col_names = T)


# Scenario 2: Colorado/UT -------------------------------------------------
# Starting and final optimal phenos set, allow rest of the params to vary. 
late_UT_CO <- NULL
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
  
  if (fitness_width == 0.65) {
    next
  }
  
  if (init_opt == 0.13 & final_opt == 0.876) {
    # Import CSV and add parameters
    one_sim_res <- read.csv(file) %>%
      mutate(sim_gens = sim_gens,
             init_corin = init_Corin,
             init_ednrb = init_EDNRB,
             init_opt = init_opt,
             final_opt = final_opt,
             fitness_width = fitness_width,
             lambda = lambda,
             replicate = replicate)
    late_UT_CO <- rbind(late_UT_CO, one_sim_res)
  }
}

write_csv(late_UT_CO, path = "results/tmp/slim_summaries/utah_and_colo.csv", col_names = T)


# Scenario 3 --------------------------------------------------------------
change_in_opt <- NULL
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
  
  if (fitness_width == 0.65) {
    next
  }
  
  if (lambda == 4 & init_opt == 0 & init_Corin == 0.05 & init_EDNRB == 0.05) {
      
    # Import CSV and add parameters
    one_sim_res <- read.csv(file) %>%
      mutate(sim_gens = sim_gens,
             init_corin = init_Corin,
             init_ednrb = init_EDNRB,
             init_opt = init_opt,
             final_opt = final_opt,
             fitness_width = fitness_width,
             lambda = lambda,
             replicate = replicate)
    change_in_opt <- rbind(change_in_opt, one_sim_res)
  }
}



write_csv(change_in_opt, path = "results/tmp/slim_summaries/changing_opt.csv", col_names = T)



# Scenario 4: UT/CO, receesive architecture -------------------------------
late_recessive_UT_CO <- NULL
for (file in list.files("slim_results/recessive_constantK/", pattern = "late.csv", full.names = T)) {
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
  
  if (fitness_width == 0.65) {
    next
  }
  
  if (init_opt == 0.13 & final_opt == 0.876) {
    # Import CSV and add parameters
    one_sim_res <- read.csv(file) %>%
      mutate(sim_gens = sim_gens,
             init_corin = init_Corin,
             init_ednrb = init_EDNRB,
             init_opt = init_opt,
             final_opt = final_opt,
             fitness_width = fitness_width,
             lambda = lambda,
             replicate = replicate)
    late_recessive_UT_CO <- rbind(late_recessive_UT_CO, one_sim_res)
  }
}

write_csv(late_recessive_UT_CO, path = "results/tmp/slim_summaries/utah_and_colo_recessive.csv", col_names = T)



# Scenario 5: UT/CO, varying K with additive genetic architecture ---------
vary_K_test <- NULL
for (file in list.files("slim_results/additive_varyK/", pattern = "late.csv", full.names = T)) {
  # parse filename for sim parameters
  elements <- str_split(basename(file), pattern = "_")[[1]]
  sim_gens <- as.numeric(str_remove(elements[1], "gens"))
  min_K <- as.numeric(str_remove(elements[2], "minK"))
  max_K <- as.numeric(str_remove(elements[3], "maxK"))
  period <- as.numeric(str_remove(elements[4], "period"))
  start_K <- as.numeric(str_remove(elements[5], "startK"))
  init_dec <- str_remove(elements[6], "initDec")
  init_Corin <- as.numeric(str_remove(elements[7], "iCorin"))
  init_EDNRB <- as.numeric(str_remove(elements[8], "iEDNRB"))
  init_opt <- as.numeric(str_remove(elements[9], "initOpt"))
  final_opt <- as.numeric(str_remove(elements[10], "finalOpt"))
  fitness_width <- as.numeric(str_remove(elements[11], "sel"))
  lambda  <- as.numeric(str_remove(elements[13], "lambda"))
  replicate <- as.numeric(str_remove(elements[14], "rep"))
  
  # if (init_opt == 0.13 & final_opt == 0.876) {
  # Import CSV and add parameters
  one_sim_res <- read.csv(file) %>%
    mutate(sim_gens = sim_gens,
           min_K = min_K,
           max_K = max_K,
           period = period,
           init_K = start_K,
           init_dec = init_dec,
           init_corin = init_Corin,
           init_ednrb = init_EDNRB,
           init_opt = init_opt,
           final_opt = final_opt,
           fitness_width = fitness_width,
           lambda = lambda,
           replicate = replicate)
  vary_K_test <- rbind(vary_K_test, one_sim_res)
  # }
}

write_csv(vary_K_test, path = "results/tmp/slim_summaries/varying_K_initial.csv", col_names = T)







