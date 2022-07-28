library(tidyverse)
library(ggrepel)
library(readxl)
library(cowplot)
library(ggokabeito)


# First, take a look at average coverage and even-ness of the final data. 
coverage <- read_delim("results/multiqc/multiqc_report_dedup_realigned_sorted_bams_data/multiqc_general_stats.txt", delim = "\t") %>% 
  mutate(sample = str_remove(Sample, ".realigned")) %>% 
  select(-Sample) %>% 
  rename_with(.cols = everything(), ~str_remove(string = .x, pattern = "QualiMap_mqc-generalstats-qualimap-")) %>% 
  pivot_longer(avg_gc:general_error_rate, names_to = "metric") %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(names_from = "metric", values_from = "value")

# Mean coverage
mean(coverage$mean_coverage) # 2.36. Not bad! Better than expected, really

# CV
sd(coverage$mean_coverage)/mean(coverage$mean_coverage) # 0.577, not so even. But perhaps not unexpected.


# "Failed" samples
# Can look for ones with coverage < 1X?
coverage %>% 
  filter(mean_coverage < 1)
# 6 samples

# Another possible metric: insert size
# The overall average is 315, which is around what we were expecting
mean(coverage$median_insert_size)
coverage %>% 
  filter(median_insert_size < 100)
# But, a few of the old museum specimens have insert sizes < 100. Might be hard to use...



# Read in the data
seqed_samples <- read_delim(file = "data/all_wtjr_sample_bamlist.txt", delim = "\t", col_names = "bam_file") %>% 
  mutate(sample = basename(bam_file)) %>% 
  mutate(sample = str_remove(sample, ".realigned.sorted.bam")) %>% 
  mutate(sample = str_remove(sample, "_realigned_sorted.bam"))


# My data
tim <- read_csv("data/labwork/possible_samples_to_extracts.csv") %>% 
  select(sample = ID_museum, year, latitude, longitude, winter_pheno, state, source = source.x) %>% 
  mutate(museum = year < 2008) %>% 
  mutate(museum = ifelse(is.na(museum), T, museum)) %>% 
  mutate(batch = "tim")

# Fix the LTW_TWO_4140 sample
tim$sample[tim$sample == "Mafalda_LTW_WYO_4140"] <- "LTW_COL_4140"
tim$year[tim$sample == "LTW_COL_4140"] <- NA
tim$latitude[tim$sample == "LTW_COL_4140"] <- 39.04814
tim$longitude[tim$sample == "LTW_COL_4140"] <- -105.56139
tim$winter_pheno[tim$sample == "LTW_COL_4140"] <- NA
tim$state[tim$sample == "LTW_COL_4140"] <- "CO"
tim$source[tim$sample == "LTW_COL_4140"] <- "LTW_COL_4140"


# Mafalda's data
mafalda <- read_excel("data/mf_thesis_table_S3_3.xlsx") %>% 
  filter(Species == "Lepus townsendii") %>%
  mutate(museum = ifelse(Origin == "Field", F, T),
         year = NA, 
         winter_pheno = NA,
         sample = str_replace(`Sample ID`, "_", "-"),
         batch = "mafalda") %>% 
  separate(`Location (Population)`, into = c("country", "state"), sep = "/") %>% 
  select(sample, year, latitude = `Lat.`, longitude = `Long.`,  winter_pheno, source = Origin, state, museum, batch) %>% 
  bind_rows(., read_csv("data/mafalda_UCM_samples.csv"))


samples <- seqed_samples %>% 
  left_join(bind_rows(tim, mafalda)) 




# Import the covariance matrix
cov_file <- "results/pcangsd/all_wtjr_samples/pca_genome_wide_GL_all_samples.cov"
C <- as.matrix(read.table(cov_file))
# Eigenvalues and eigenvectors
e <- eigen(C)

e$values/sum(e$values)

# Get just the eigenvectors, which are apparently our principal components
# (but not our PC scores?)
# This doesn't really make sense to me: I though I needed to project the original genotypes
# Back into the PC space?

pcs <- as_tibble(e$vectors) %>% 
  rename_with(str_replace, pattern = "V", replacement = "PC") 

pca_res <- bind_cols(samples, pcs) %>% 
  mutate(state = ifelse(state == "Colorado", "CO", state)) %>%
  mutate(state = ifelse(state == "utah", "UT", state))
# Some prelim plots

# PCA: color by altitude,
a <- ggplot(aes(x = PC1, y = PC2, color = batch, label = sample), data = pca_res) +
  geom_label_repel(max.overlaps = 30, min.segment.length = 0.2, force = 4, force_pull = 0.5) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_discrete(type = c("orange", "dodgerblue3")) +
  theme_cowplot() # +
  #ggtitle("PCAngsd, genome-wide, color by elevation")

a

b <- ggplot(aes(x = PC1, y = PC2, color = museum, label = sample), data = pca_res) +
  #geom_label_repel(max.overlaps = 30, min.segment.length = 0.2, force = 4, force_pull = 0.5) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_discrete(type = c("orange", "dodgerblue3")) +
  theme_cowplot()

# +
#ggtitle("PCAngsd, genome-wide, color by elevation")



cor.test(pca_res$PC1, pca_res$latitude)

ggplot(aes(x = PC1, y = PC2, color = state, label = sample), data = pca_res) +
  #geom_label_repel(max.overlaps = 30, min.segment.length = 0.2, force = 4, force_pull = 0.5) +
  geom_point(size = 7, alpha = 0.8) +
  theme_cowplot() +
  scale_color_okabe_ito() +
  xlab("PC1- 8.98%") +
  ylab("PC2- 1.32%") +
  theme(legend.position = c(0.4, 0.06),
        legend.direction = "horizontal",
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.title = element_text(size = 26),
        axis.text = element_text(size = 22),
        axis.line = element_line(size = 1.5))
ggsave(filename = "results/genomic_pca_state_inset.pdf")


ggplot(aes(x = PC1, y = latitude), data = pca_res) +
  geom_point(aes(color = state), size = 4, alpha = 0.8) +
  #geom_smooth(method = "lm") +
  theme_cowplot() +
  scale_color_okabe_ito() +
  xlab("PC1- 8.98%") +
  ylab("Latitude")
ggsave(filename = "results/pc1_by_latitude.pdf", height = 6, width = 8, units = "in")










