
# CO/UT- additive, constant K ---------------------------------------------
additive_constant_K <- NULL
for (file in list.files(
  "slim_results_final/additive_constantK/",
  pattern = "late.csv",
  full.names = T
)) {
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
    mutate(
      sim_gens = sim_gens,
      init_corin = init_Corin,
      init_ednrb = init_EDNRB,
      init_opt = init_opt,
      final_opt = final_opt,
      fitness_width = fitness_width,
      lambda = lambda,
      replicate = replicate
    )
  additive_constant_K <- rbind(additive_constant_K, one_sim_res)
}


# Check for proper number of replicates for each combination of parameters
additive_constant_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>% 
  filter(n != 30)

# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 4*4*3*2 == 96
additive_constant_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>%
  ungroup() %>% 
  dim(.) %>% 
  .[1] == 96

# Proper number in each category
length(unique(additive_constant_K$init_corin)) == 4
length(unique(additive_constant_K$init_ednrb)) == 4
length(unique(additive_constant_K$init_opt)) == 1
length(unique(additive_constant_K$final_opt)) == 1
length(unique(additive_constant_K$sim_gens)) == 1
length(unique(additive_constant_K$fitness_width)) == 3
length(unique(additive_constant_K$lambda)) == 2
length(unique(additive_constant_K$replicate)) == 30


# Save results
write_csv(additive_constant_K, path = "results/slim_summaries/additive_constantK.csv", col_names = T)


# UT/CO, recessive constant K -----------------------------------------------------
recessive_constant_K <- NULL
for (file in list.files(
  "slim_results_final/recessive_constantK/",
  pattern = "late.csv",
  full.names = T
)) {
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
    mutate(
      sim_gens = sim_gens,
      init_corin = init_Corin,
      init_ednrb = init_EDNRB,
      init_opt = init_opt,
      final_opt = final_opt,
      fitness_width = fitness_width,
      lambda = lambda,
      replicate = replicate
    )
  recessive_constant_K <- rbind(recessive_constant_K, one_sim_res)
}



# Check for proper number of replicates for each combination of parameters
recessive_constant_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>% 
  filter(n != 30)

# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 4*4*3*2 == 96
recessive_constant_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>%
  ungroup() %>% 
  dim(.) %>% 
  .[1] == 96

# Proper number in each category
length(unique(recessive_constant_K$init_corin)) == 4
length(unique(recessive_constant_K$init_ednrb)) == 4
length(unique(recessive_constant_K$init_opt)) == 1
length(unique(recessive_constant_K$final_opt)) == 1
length(unique(recessive_constant_K$sim_gens)) == 1
length(unique(recessive_constant_K$fitness_width)) == 3
length(unique(recessive_constant_K$lambda)) == 2
length(unique(recessive_constant_K$replicate)) == 30


# Save results
write_csv(recessive_constant_K, path = "results/slim_summaries/recessive_constantK.csv", col_names = T)


# UT/CO, additive, varying K ----------------------------------------------
additive_vary_K <- NULL
for (file in list.files("slim_results_final/additive_varyK/", pattern = "late.csv", full.names = T)) {
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
  additive_vary_K <- rbind(additive_vary_K, one_sim_res)

}


# Check for proper number of replicates for each combination of parameters
additive_vary_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda, init_dec) %>% 
  tally() %>% 
  filter(n != 30)

# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 2 initial directions
# 4*4*3*2*2 == 192
additive_vary_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda, init_dec) %>% 
  tally() %>%
  ungroup() %>% 
  dim(.) %>% 
  .[1] == 192


# Proper number in each category
length(unique(additive_vary_K$init_corin)) == 4
length(unique(additive_vary_K$init_ednrb)) == 4
length(unique(additive_vary_K$init_opt)) == 1
length(unique(additive_vary_K$final_opt)) == 1
length(unique(additive_vary_K$sim_gens)) == 1
length(unique(additive_vary_K$fitness_width)) == 3
length(unique(additive_vary_K$lambda)) == 2
length(unique(additive_vary_K$replicate)) == 30
length(unique(additive_vary_K$init_dec)) == 2


# Save results
write_csv(additive_vary_K, path = "results/slim_summaries/additive_varyK.csv", col_names = T)



