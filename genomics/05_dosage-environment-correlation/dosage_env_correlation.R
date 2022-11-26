

library(corrplot)
library(raster)
library(tidyverse)



# Load sample information
wtjr <- read_csv("raw_data/genetics/samples_info.csv")


# Load environmental data
env_variables <- raster::stack("processed_data/pheno_predictors_millsetal2018.tif")
names(env_variables) <- c("snow.cover", "bio_2", "bio_3")


# Extract environmental data and combine
env_preds <- raster::extract(env_variables, 
                             y = sp::SpatialPoints(data.frame(x = as.numeric(wtjr$longitude),
                                                              y = as.numeric(wtjr$latitude)), 
                                                   proj4string = CRS(proj4string(env_variables)))) %>% 
  as.data.frame(.)

wtjr_env <- wtjr %>% 
  bind_cols(env_preds)


# Read in genotype dosages of white allele
read_in_white_dosage <- function(filename) {
  white_dosage <- read_csv(filename) %>% 
    dplyr::select(-scaffold, -position) %>% 
    t(.) %>% 
    as_tibble(., rownames = "sample") 
  colnames(white_dosage) <- white_dosage[1,]
  colnames(white_dosage)[1] <- "sample"
  white_dosage <- white_dosage[-1,]
  
  # Return
  white_dosage %>% 
    mutate(across(-sample, as.numeric))
}

# All associated SNPs for which we have dosages,
# in terms of white dosage
corin_white_dosage <- read_in_white_dosage("genomics/results/pcangsd/95ind_10X_filter/corin_dosages_white.csv")
asip_white_dosage <- read_in_white_dosage("genomics/results/pcangsd/95ind_10X_filter/asip_dosages_white.csv")
ednrb_white_dosage <- read_in_white_dosage("genomics/results/pcangsd/95ind_10X_filter/ednrb_dosages_white.csv")
scaff380_white_dosage <- read_in_white_dosage("genomics/results/pcangsd/95ind_10X_filter/scaff380_dosages_white.csv")


# Get top associated SNP across each site
corin_top_white <- corin_white_dosage %>% 
  dplyr::select(sample, top = `342_46966737`) %>% 
  mutate(gene = "corin")

asip_top_white <- asip_white_dosage %>% 
  dplyr::select(sample, top = `245_24236852`) %>% 
  mutate(gene = "asip")

ednrb_top_white <- ednrb_white_dosage %>% 
  dplyr::select(sample, top = `311_3413288`) %>% 
  mutate(gene = "ednrb")

scaff380_top_white <- scaff380_white_dosage %>% 
  dplyr::select(sample, top = `380_35303819`) %>% 
  mutate(gene = "scaff380")

dosage_top <- bind_rows(asip_top_white, corin_top_white, ednrb_top_white, scaff380_top_white)

# Combine it all
wtjr_white_dosages <- wtjr_env %>% 
  left_join(dosage_top) 


for_top_cor <- wtjr_white_dosages %>% 
  dplyr::select(sample, snow.cover, bio_2, bio_3, gene, top) %>% 
  pivot_wider(names_from = gene, values_from = top) %>% 
  dplyr::select(-sample) %>% 
  relocate(asip:scaff380)


# Make corrplot for supplemental material. 
png("genomics/results/figures/supplemental/GEA_corrplot_95ind.png", width = 800, height = 800)
corr.mat <- cor(for_top_cor)
pmat <- cor.mtest(for_top_cor, conf.level = 0.95)
corr_plot <- corrplot(corr.mat, p.mat = pmat$p, sig.level = 0.05,
                      diag = F, method = 'circle', insig='blank', title = " ",
                      tl.col = "black",
                      col = COL2("PRGn"),
                      type = "lower")$corrPos
text(corr_plot$x, corr_plot$y, round(corr_plot$corr, 2), )
dev.off()


