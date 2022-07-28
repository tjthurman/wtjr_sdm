library(cowplot)
library(IRanges)
library(biomaRt)
library(tidyverse)


# Import data
res_folder <- "results/angsd_Fst/colorado_samples_bamlist_NDK_samples_bamlist/wtjr_COL-NDK_20220503/fst/winSize50000/winStep50000"

# Function that loads and compiles Fsr data from a given folder.
# Pretty untested, relies on a very specific naming scheme of files to work (the scheme I am using). 
compile_angsd_fst_results <- function(res_folder) {
  fsts <- tibble(NULL)
  window_size <- as.numeric(str_remove(str_split(res_folder, "/")[[1]][str_detect(str_split(res_folder, "/")[[1]], "winSize")], "winSize"))
  step_size <- as.numeric(str_remove(str_split(res_folder, "/")[[1]][str_detect(str_split(res_folder, "/")[[1]], "winStep")], "winStep"))
  
  for (file in list.files(res_folder, pattern = "fst_")) {
    scaffold <- str_remove(str_split(file, pattern = "_")[[1]][2], ".txt")
    
    fst  <- read_delim(file.path(res_folder, file), delim = "\t", skip = 1, 
                       col_names = c("region", "scaff", "midpoint", "sites", "fst"), 
                       col_types = "cciin" ) %>% 
      mutate(scaffold = scaffold) %>% 
      dplyr::select(scaffold, scaff, region, sites, midpoint, fst)
    
    if (sum(fst$scaffold == fst$scaff) != dim(fst)[1]) {
      stop("scaffold names don't match")
    }
    
    fsts <- rbind(fsts, dplyr::select(fst, -scaff))
  }
  
  return(fsts)
}

# Compile results form original way, using the default options in angsd
all_fst <- compile_angsd_fst_results(res_folder)

# Get-genome-wide mean
mean(all_fst$fst)
max(all_fst$fst)
# Get genome-wide mean for all windows with at least 5k sites in the 50kb window
all_fst %>% 
  filter(sites > 500) %>% 
  summarize(mean(fst))


# ight be good to do the genome-wide estimate in angsd
# Would need to cat all the beagle gls, run fst index, and then run fst stats. Would take a long time...

# Doing a genome wide estimate may not work. trying to concatenate all the per-sacffold .saf files lead to some massive,
# massive files. Instead, ran ran realSFS stats for each scaffold to get the per-scaffold .fst values. Can then take the mean of those
# perhaps weighting by the number of observations per scaffold?

raw <- read_lines("results/angsd_Fst/colorado_samples_bamlist_NDK_samples_bamlist/wtjr_COL-NDK_20220503/fst/genome_fst_stats_res_stderr.txt")


scaffolds_1 <- raw[seq(from = 1, to = length(raw), by = 3)] %>% 
  str_extract("[:digit:]+") %>% 
  as.numeric()
scaffold_2 <- raw[seq(from = 2, to = length(raw), by = 3)] %>% 
  str_extract("[:digit:]+") %>% 
  as.numeric()
fst_info <- raw[seq(from = 3, to = length(raw), by = 3)] %>% 
  str_remove("\\\t-> ") %>% 
  str_remove("FST\\.Unweight") %>%
  str_remove(" Fst\\.Weight") %>% 
  str_remove("nObs\\:") %>% 
  str_remove_all("\\[") %>% 
  str_remove_all("\\]") %>% 
  str_remove_all(" ") 

per_scaffold_fst <- bind_cols(scaffold = scaffolds_1,
          fst_string = fst_info) %>%
  separate(col = "fst_string", 
           into = c("nObs", "unweight_fst", "weight_fst"),
           sep = ":", 
           convert = T) 


# Now, can calculate mean genome-wide Fst a few ways. Can just take the naive mean of the per-scaffold weighted and unweighted Fsts:
per_scaffold_fst %>% 
  pivot_longer(unweight_fst:weight_fst) %>% 
  group_by(name) %>% 
  summarize(mean = mean(value, na.rm = T),
           sd = sd(value, na.rm = T),
           quantile05 = quantile(value, na.rm = T, probs = 0.05),
           quantile95 = quantile(value, na.rm = T, probs = 0.95))


per_scaffold_fst %>% 
  pivot_longer(unweight_fst:weight_fst) %>% 
  group_by(name) %>% 
  summarize(weight_mean = weighted.mean(value, na.rm = T, w = nObs),
            sd = sd(value, na.rm = T),
            quantile05 = quantile(value, na.rm = T, probs = 0.05),
            quantile95 = quantile(value, na.rm = T, probs = 0.95))


# A final way: results of SFS print

print_res <- read_delim("results/angsd_Fst/colorado_samples_bamlist_NDK_samples_bamlist/wtjr_COL-NDK_20220503/fst/genome_fst_print.txt",
                        delim  = "\t", col_names = c("scaffold", "position", "a", "b"))

x <- print_res %>% 
  mutate(sum =  a + b) %>% 
  mutate(ratio = a/(sum)) %>% 
  summarize(weighted_fst = sum(a, na.rm = T)/sum(sum, na.rm = T),
            weighted_fst_2 = sum(a, na.rm = T)/(sum(a, na.rm = T) + sum(b, na.rm = T)),
            unweighted_fst = mean(ratio, na.rm = T))
