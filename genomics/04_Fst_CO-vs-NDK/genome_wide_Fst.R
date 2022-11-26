library(tidyverse)

# A final way: results of SFS print

# This file is generated from the `print_Fst_for_genome_wide.sh` script and is too large to store on Git. 
print_res <- read_delim("genomics/results/angsd_Fst/colorado_samples_bamlist_NDK_samples_bamlist/wtjr_COL-NDK_20221113/fst/genome_fst_print.txt",
                        delim  = "\t", col_names = c("scaffold", "position", "a", "b"))

fst <- print_res %>% 
  mutate(sum =  a + b) %>% 
  mutate(ratio = a/(sum)) %>% 
  summarize(weighted_fst = sum(a, na.rm = T)/sum(sum, na.rm = T))
