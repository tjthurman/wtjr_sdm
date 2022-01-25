library(tidyverse)
library(assertthat)


# CO/UT- additive, constant K ---------------------------------------------
additive_constant_K <- NULL
for (file in list.files(
  "slim_results_final_80gens_consVE/additive_constantK/",
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
not30 <- additive_constant_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>% 
  filter(n != 30)

assert_that(dim(not30)[1] == 0, msg = "not all replicates present for additive_constant_K")

# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 4*4*3*2 == 96
param_combos <- additive_constant_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>%
  ungroup() 

assert_that(dim(param_combos)[1] == 96, msg = "Incorrect number of parameter combinations for additive_constant_K")


# Proper number in each category
assert_that(length(unique(additive_constant_K$init_corin)) == 4, msg = "Incorrect initial CORIN for additive_constant_K")
assert_that(length(unique(additive_constant_K$init_ednrb)) == 4, msg = "Incorrect initial EDNRB for additive_constant_K")
assert_that(length(unique(additive_constant_K$init_opt)) == 1, msg = "Incorrect initial pheno opt for additive_constant_K")
assert_that(length(unique(additive_constant_K$final_opt)) == 1, msg = "Incorrect final pheno opt for additive_constant_K")
assert_that(length(unique(additive_constant_K$sim_gens)) == 1, msg = "Incorrect generations for additive_constant_K")
assert_that(length(unique(additive_constant_K$fitness_width)) == 3, msg = "Incorrect fitness widths for additive_constant_K")
assert_that(length(unique(additive_constant_K$lambda)) == 2, msg = "Incorrect lambdasfor additive_constant_K")
assert_that(length(unique(additive_constant_K$replicate)) == 30, msg = "Incorrect replicates for additive_constant_K")

# Save results
write_csv(additive_constant_K, path = "results/slim_summaries_80gens_consVE/additive_constantK.csv", col_names = T)

# rm intermediate files
rm(not30, param_combos)

# UT/CO, recessive constant K -----------------------------------------------------
recessive_constant_K <- NULL
for (file in list.files(
  "slim_results_final_80gens/recessive_constantK/",
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
not30 <- recessive_constant_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>% 
  filter(n != 30)

assert_that(dim(not30)[1] == 0, msg = "not all replicates present for recessive_constant_K")

# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 4*4*3*2 == 96
param_combos <- recessive_constant_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>%
  ungroup() 

assert_that(dim(param_combos)[1] == 96, msg = "Incorrect number of parameter combinations for recessive_constant_K")


# Proper number in each category
assert_that(length(unique(recessive_constant_K$init_corin)) == 4, msg = "Incorrect initial CORIN for recessive_constant_K")
assert_that(length(unique(recessive_constant_K$init_ednrb)) == 4, msg = "Incorrect initial EDNRB for recessive_constant_K")
assert_that(length(unique(recessive_constant_K$init_opt)) == 1, msg = "Incorrect initial pheno opt for recessive_constant_K")
assert_that(length(unique(recessive_constant_K$final_opt)) == 1, msg = "Incorrect final pheno opt for recessive_constant_K")
assert_that(length(unique(recessive_constant_K$sim_gens)) == 1, msg = "Incorrect generations for recessive_constant_K")
assert_that(length(unique(recessive_constant_K$fitness_width)) == 3, msg = "Incorrect fitness widths for recessive_constant_K")
assert_that(length(unique(recessive_constant_K$lambda)) == 2, msg = "Incorrect lambdasfor recessive_constant_K")
assert_that(length(unique(recessive_constant_K$replicate)) == 30, msg = "Incorrect replicates for recessive_constant_K")

# Save results
write_csv(recessive_constant_K, path = "results/slim_summaries_80gens/recessive_constantK.csv", col_names = T)

# rm intermediate files
rm(not30, param_combos)


# UT/CO, additive, varying K ----------------------------------------------
additive_vary_K <- NULL
for (file in list.files("slim_results_final_80gens/additive_varyK/", pattern = "late.csv", full.names = T)) {
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
not30 <- additive_vary_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda, init_dec) %>% 
  tally() %>% 
  filter(n != 30)

assert_that(dim(not30)[1] == 0, msg = "not all replicates present for additive_vary_K")



# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 2 initial directions
# 4*4*3*2*2 == 192
param_combos <- additive_vary_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda, init_dec) %>% 
  tally() %>%
  ungroup() 

assert_that(dim(param_combos)[1] == 192, msg = "Incorrect number of parameter combinations for additive_vary_K")


# Proper number in each category
assert_that(length(unique(additive_vary_K$init_corin)) == 4, msg = "Incorrect initial CORIN for additive_vary_K")
assert_that(length(unique(additive_vary_K$init_ednrb)) == 4, msg = "Incorrect initial EDNRB for additive_vary_K")
assert_that(length(unique(additive_vary_K$init_opt)) == 1, msg = "Incorrect initial pheno opt for additive_vary_K")
assert_that(length(unique(additive_vary_K$final_opt)) == 1, msg = "Incorrect final pheno opt for additive_vary_K")
assert_that(length(unique(additive_vary_K$sim_gens)) == 1, msg = "Incorrect generations for additive_vary_K")
assert_that(length(unique(additive_vary_K$fitness_width)) == 3, msg = "Incorrect fitness widths for additive_vary_K")
assert_that(length(unique(additive_vary_K$lambda)) == 2, msg = "Incorrect lambdasfor additive_vary_K")
assert_that(length(unique(additive_vary_K$replicate)) == 30, msg = "Incorrect replicates for additive_vary_K")
assert_that(length(unique(additive_vary_K$init_dec)) == 2, msg = "Incorrect replicates for additive_vary_K")

# Save results
write_csv(additive_vary_K, path = "results/slim_summaries_80gens/additive_varyK.csv", col_names = T)

# rm intermediate files
rm(not30, param_combos)

# Recessive, 1 locus (SSH) ------------------------------------------------
SSH_constant_K <- NULL
for (file in list.files(
  "slim_results_final_80gens/recessive_constantK_SSH//",
  pattern = "late.csv",
  full.names = T
)) {
  # parse filename for sim parameters
  elements <- str_split(basename(file), pattern = "_")[[1]]
  sim_gens <- as.numeric(str_remove(elements[1], "gens"))
  init_Corin <- as.numeric(str_remove(elements[4], "iCorin"))
  init_opt <- as.numeric(str_remove(elements[5], "initOpt"))
  final_opt <- as.numeric(str_remove(elements[6], "finalOpt"))
  fitness_width <- as.numeric(str_remove(elements[7], "sel"))
  lambda  <- as.numeric(str_remove(elements[9], "lambda"))
  replicate <- as.numeric(str_remove(elements[10], "rep"))
  
  # Import CSV and add parameters
  one_sim_res <- read.csv(file) %>%
    mutate(
      sim_gens = sim_gens,
      init_corin = init_Corin,
      init_opt = init_opt,
      final_opt = final_opt,
      fitness_width = fitness_width,
      lambda = lambda,
      replicate = replicate
    )
  SSH_constant_K <- rbind(SSH_constant_K, one_sim_res)
}



# Check for proper number of replicates for each combination of parameters
# Check for proper number of replicates for each combination of parameters
not30 <- SSH_constant_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>% 
  filter(n != 30)

assert_that(dim(not30)[1] == 0, msg = "not all replicates present for SSH_constant_K")

# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 4*3*2 == 24
param_combos <- SSH_constant_K %>% 
  filter(generation == 1) %>% 
  group_by(init_corin,, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>%
  ungroup()

assert_that(dim(param_combos)[1] == 24, msg = "Incorrect number of parameter combinations for SSH_constant_K")


# Proper number in each category
# Proper number in each category
assert_that(length(unique(SSH_constant_K$init_corin)) == 4, msg = "Incorrect initial CORIN for SSH_constant_K")
assert_that(length(unique(SSH_constant_K$init_opt)) == 1, msg = "Incorrect initial pheno opt for SSH_constant_K")
assert_that(length(unique(SSH_constant_K$final_opt)) == 1, msg = "Incorrect final pheno opt for SSH_constant_K")
assert_that(length(unique(SSH_constant_K$sim_gens)) == 1, msg = "Incorrect generations for SSH_constant_K")
assert_that(length(unique(SSH_constant_K$fitness_width)) == 3, msg = "Incorrect fitness widths for SSH_constant_K")
assert_that(length(unique(SSH_constant_K$lambda)) == 2, msg = "Incorrect lambdasfor SSH_constant_K")
assert_that(length(unique(SSH_constant_K$replicate)) == 30, msg = "Incorrect replicates for SSH_constant_K")

# Save results
write_csv(SSH_constant_K, path = "results/slim_summaries_80gens/SSH_constantK.csv", col_names = T)

# rm intermediate files
rm(not30, param_combos)

