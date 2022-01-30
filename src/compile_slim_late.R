####################
# Project: WTJR SDMs
# Author: Timothy Thurman
# Purpose: Compile individual SLiM late files
####################


# Load packages -----------------------------------------------------------
library(tidyverse)
library(assertthat)


# List folders ------------------------------------------------------------
additive_constantK_2locus_folder <- "results/slim_ind_sims/additive_constantK_2locus/"
recessive_constantK_2locus_folder <- "results/slim_ind_sims/recessive_constantK_2locus/"
recessive_constantK_1locus_folder <- "results/slim_ind_sims/recessive_constantK_1locus/"
additive_varyK_2locus_folder <- "results/slim_ind_sims/additive_varyK_2locus/"


# additive, constant K, 2 locus, files ---------------------------------------------
additive_consK_2locus_late <- NULL
for (file in list.files(
  additive_constantK_2locus_folder,
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
  additive_consK_2locus_late <- rbind(additive_consK_2locus_late, one_sim_res)
  rm(one_sim_res, elements, file, final_opt, fitness_width, init_Corin, init_EDNRB, init_opt, lambda, replicate, sim_gens)
}


# Check for proper number of replicates for each combination of parameters
not30 <- additive_consK_2locus_late %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>% 
  filter(n != 30)

assert_that(dim(not30)[1] == 0, msg = "not all replicates present for additive_consK_2locus_late")

# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 4*4*3*2 == 96
param_combos <- additive_consK_2locus_late %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>%
  ungroup() 

assert_that(dim(param_combos)[1] == 96, msg = "Incorrect number of parameter combinations for additive_consK_2locus_late")


# Proper number in each category
assert_that(length(unique(additive_consK_2locus_late$init_corin)) == 4, msg = "Incorrect initial CORIN for additive_consK_2locus_late")
assert_that(length(unique(additive_consK_2locus_late$init_ednrb)) == 4, msg = "Incorrect initial EDNRB for additive_consK_2locus_late")
assert_that(length(unique(additive_consK_2locus_late$init_opt)) == 1, msg = "Incorrect initial pheno opt for additive_consK_2locus_late")
assert_that(length(unique(additive_consK_2locus_late$final_opt)) == 1, msg = "Incorrect final pheno opt for additive_consK_2locus_late")
assert_that(length(unique(additive_consK_2locus_late$sim_gens)) == 1, msg = "Incorrect generations for additive_consK_2locus_late")
assert_that(length(unique(additive_consK_2locus_late$fitness_width)) == 3, msg = "Incorrect fitness widths for additive_consK_2locus_late")
assert_that(length(unique(additive_consK_2locus_late$lambda)) == 2, msg = "Incorrect lambdasfor additive_consK_2locus_late")
assert_that(length(unique(additive_consK_2locus_late$replicate)) == 30, msg = "Incorrect replicates for additive_consK_2locus_late")

# Save results
write_csv(additive_consK_2locus_late, path = "results/slim_summaries/additive_constantK_2locus_late.csv", col_names = T)

# rm intermediate files
rm(not30, param_combos)


# recessive, constant K, 2 locus, late files ---------------------------------------------
recessive_consK_2locus_late <- NULL
for (file in list.files(
  recessive_constantK_2locus_folder,
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
  recessive_consK_2locus_late <- rbind(recessive_consK_2locus_late, one_sim_res)
  rm(one_sim_res, elements, sim_gens, init_Corin, init_EDNRB, init_opt, final_opt, fitness_width, lambda, replicate)
}


# Check for proper number of replicates for each combination of parameters
not30 <- recessive_consK_2locus_late %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>% 
  filter(n != 30)

assert_that(dim(not30)[1] == 0, msg = "not all replicates present for recessive_consK_2locus_late")

# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 4*4*3*2 == 96
param_combos <- recessive_consK_2locus_late %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>%
  ungroup() 

assert_that(dim(param_combos)[1] == 96, msg = "Incorrect number of parameter combinations for recessive_consK_2locus_late")


# Proper number in each category
assert_that(length(unique(recessive_consK_2locus_late$init_corin)) == 4, msg = "Incorrect initial CORIN for recessive_consK_2locus_late")
assert_that(length(unique(recessive_consK_2locus_late$init_ednrb)) == 4, msg = "Incorrect initial EDNRB for recessive_consK_2locus_late")
assert_that(length(unique(recessive_consK_2locus_late$init_opt)) == 1, msg = "Incorrect initial pheno opt for recessive_consK_2locus_late")
assert_that(length(unique(recessive_consK_2locus_late$final_opt)) == 1, msg = "Incorrect final pheno opt for recessive_consK_2locus_late")
assert_that(length(unique(recessive_consK_2locus_late$sim_gens)) == 1, msg = "Incorrect generations for recessive_consK_2locus_late")
assert_that(length(unique(recessive_consK_2locus_late$fitness_width)) == 3, msg = "Incorrect fitness widths for recessive_consK_2locus_late")
assert_that(length(unique(recessive_consK_2locus_late$lambda)) == 2, msg = "Incorrect lambdasfor recessive_consK_2locus_late")
assert_that(length(unique(recessive_consK_2locus_late$replicate)) == 30, msg = "Incorrect replicates for recessive_consK_2locus_late")

# Save results
write_csv(recessive_consK_2locus_late, path = "results/slim_summaries/recessive_constantK_2locus_late.csv", col_names = T)

# rm intermediate files
rm(not30, param_combos)


# recessive, constant K, 1 locus, late files ---------------------------------------------
recessive_consK_1locus_late <- NULL
for (file in list.files(
  recessive_constantK_1locus_folder,
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
  recessive_consK_1locus_late <- rbind(recessive_consK_1locus_late, one_sim_res)
  rm(one_sim_res, elements, sim_gens, init_Corin, init_opt, final_opt, fitness_width, lambda, replicate)
}



# Check for proper number of replicates for each combination of parameters
not30 <- recessive_consK_1locus_late %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>% 
  filter(n != 30)

assert_that(dim(not30)[1] == 0, msg = "not all replicates present for recessive_consK_1locus_late")

# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 4*3*2 == 24
param_combos <- recessive_consK_1locus_late %>% 
  filter(generation == 1) %>% 
  group_by(init_corin,, init_opt, final_opt, fitness_width, lambda) %>% 
  tally() %>%
  ungroup()

assert_that(dim(param_combos)[1] == 24, msg = "Incorrect number of parameter combinations for recessive_consK_1locus_late")


# Proper number in each category
# Proper number in each category
assert_that(length(unique(recessive_consK_1locus_late$init_corin)) == 4, msg = "Incorrect initial CORIN for recessive_consK_1locus_late")
assert_that(length(unique(recessive_consK_1locus_late$init_opt)) == 1, msg = "Incorrect initial pheno opt for recessive_consK_1locus_late")
assert_that(length(unique(recessive_consK_1locus_late$final_opt)) == 1, msg = "Incorrect final pheno opt for recessive_consK_1locus_late")
assert_that(length(unique(recessive_consK_1locus_late$sim_gens)) == 1, msg = "Incorrect generations for recessive_consK_1locus_late")
assert_that(length(unique(recessive_consK_1locus_late$fitness_width)) == 3, msg = "Incorrect fitness widths for recessive_consK_1locus_late")
assert_that(length(unique(recessive_consK_1locus_late$lambda)) == 2, msg = "Incorrect lambdasfor recessive_consK_1locus_late")
assert_that(length(unique(recessive_consK_1locus_late$replicate)) == 30, msg = "Incorrect replicates for recessive_consK_1locus_late")

# Save results
write_csv(recessive_consK_1locus_late, path = "results/slim_summaries/recessive_constantK_1locus_late.csv", col_names = T)

# rm intermediate files
rm(not30, param_combos)

# additive, vary K, 2 locus, late files ---------------------------------------------
additive_varyK_2locus_late <- NULL
for (file in list.files(additive_varyK_2locus_folder, 
                        pattern = "late.csv", 
                        full.names = T)) {
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
  additive_varyK_2locus_late <- rbind(additive_varyK_2locus_late, one_sim_res)
  rm(one_sim_res, elements, sim_gens, min_K, max_K, period, start_K, init_dec, init_Corin, init_EDNRB, init_opt, final_opt, fitness_width, lambda, replicate)
  
}


# Check for proper number of replicates for each combination of parameters
not30 <- additive_varyK_2locus_late %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda, init_dec) %>% 
  tally() %>% 
  filter(n != 30)

assert_that(dim(not30)[1] == 0, msg = "not all replicates present for additive_varyK_2locus_late")


# Check for proper number of combo of params
# 4 for each initial allele freq
# 3 selections, 2 lambdas
# 2 initial directions
# 4*4*3*2*2 == 192
param_combos <- additive_varyK_2locus_late %>% 
  filter(generation == 1) %>% 
  group_by(init_corin, init_ednrb, init_opt, final_opt, fitness_width, lambda, init_dec) %>% 
  tally() %>%
  ungroup() 

assert_that(dim(param_combos)[1] == 192, msg = "Incorrect number of parameter combinations for additive_varyK_2locus_late")


# Proper number in each category
assert_that(length(unique(additive_varyK_2locus_late$init_corin)) == 4, msg = "Incorrect initial CORIN for additive_varyK_2locus_late")
assert_that(length(unique(additive_varyK_2locus_late$init_ednrb)) == 4, msg = "Incorrect initial EDNRB for additive_varyK_2locus_late")
assert_that(length(unique(additive_varyK_2locus_late$init_opt)) == 1, msg = "Incorrect initial pheno opt for additive_varyK_2locus_late")
assert_that(length(unique(additive_varyK_2locus_late$final_opt)) == 1, msg = "Incorrect final pheno opt for additive_varyK_2locus_late")
assert_that(length(unique(additive_varyK_2locus_late$sim_gens)) == 1, msg = "Incorrect generations for additive_varyK_2locus_late")
assert_that(length(unique(additive_varyK_2locus_late$fitness_width)) == 3, msg = "Incorrect fitness widths for additive_varyK_2locus_late")
assert_that(length(unique(additive_varyK_2locus_late$lambda)) == 2, msg = "Incorrect lambdasfor additive_varyK_2locus_late")
assert_that(length(unique(additive_varyK_2locus_late$replicate)) == 30, msg = "Incorrect replicates for additive_varyK_2locus_late")
assert_that(length(unique(additive_varyK_2locus_late$init_dec)) == 2, msg = "Incorrect replicates for additive_varyK_2locus_late")

# Save results
write_csv(additive_varyK_2locus_late, path = "results/slim_summaries/additive_varyK_2locus_late.csv", col_names = T)

# rm intermediate files
rm(not30, param_combos)
