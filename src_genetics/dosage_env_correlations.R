library(readxl)
library(raster)
library(GGally)
library(tidyverse)
library(broom)
library(factoextra)
library(cowplot)

# Get samples from the bamlist
samples <- read_delim(file = "data/all_wtjr_sample_bamlist.txt", delim = "\t", col_names = "bam_file") %>% 
  mutate(sample = basename(bam_file)) %>% 
  mutate(sample = str_remove(sample, ".realigned.sorted.bam")) %>% 
  mutate(sample = str_remove(sample, "_realigned_sorted.bam"))


# My data
tim <- read_csv("data/labwork/possible_samples_to_extracts.csv") %>% 
  select(sample = ID_museum, year, latitude, longitude, winter_pheno, state, source = source.x) %>% 
  mutate(museum = year < 2008) %>% 
  mutate(museum = ifelse(is.na(museum), T, museum)) %>% 
  mutate(batch = "tim")

# Mafalda's data
mafalda <- read_excel("data/mf_thesis_table_S3_3.xlsx") %>% 
  filter(Species == "Lepus townsendii") %>%
  mutate(museum = ifelse(Origin == "Field", F, T),
         year = NA, 
         winter_pheno = NA,
         sample = str_replace(`Sample ID`, "_", "-"),
         batch = "mafalda") %>% 
  separate(`Location (Population)`, into = c("country", "state"), sep = "/") %>% 
  select(sample, year, latitude = `Lat.`, longitude = `Long.`,  winter_pheno, source = Origin, museum, batch) %>% 
  bind_rows(., read_csv("data/mafalda_UCM_samples.csv"))


# Mafalda UCM
# Some samples that she had originally sequenced but excluded from her thesis for various reasons, I guess?



corin_dosages <- read_tsv("results/pcangsd/all_wtjr_samples/corin_dosages.tsv", col_names = "corin")
asip_dosages <- read_tsv("results/pcangsd/all_wtjr_samples/asip_dosages.tsv", col_names = "asip")
ednrb_dosages <- read_tsv("results/pcangsd/all_wtjr_samples/ednrb_dosages.tsv", col_names = "ednrb")


# Combine it all
wtjr <- samples %>% 
  left_join(bind_rows(tim, mafalda)) %>% 
  bind_cols(., corin_dosages, asip_dosages, ednrb_dosages)

wtjr_locs <- wtjr %>% 
  filter(!is.na(latitude),
         !is.na(longitude))


wtjr %>% 
  pivot_longer(corin:ednrb, names_to = "gene", values_to = "dosage") %>% 
  ggplot(aes(x = latitude, y = dosage, color = gene)) +
  geom_smooth(method = "lm") +
  geom_point() 




# Load the mills et al predictors
env_variables <- raster::stack("../wtjr_sdm/processed_data/pheno_predictors_millsetal2018.tif")
names(env_variables) <- c("snow.cover", "bio_2", "bio_3")

env_preds <- raster::extract(env_variables, y = sp::SpatialPoints(data.frame(x = as.numeric(wtjr_locs$longitude),
                                                                 y = as.numeric(wtjr_locs$latitude)), 
                                                      proj4string = CRS(proj4string(env_variables)))) %>% 
  as.data.frame(.)

z <- wtjr_locs %>% 
  bind_cols(env_preds)

cor.test(z$bio_2, z$corin)

x <- z %>% 
  dplyr::select(corin:bio_3) %>% 
  ggpairs(.)

ggsave(x, file = "results/dosage_env_cor.pdf", height = 7, width = 10, unit = "in")

z %>% 
  dplyr::select(corin:ednrb) %>% 
  ggpairs(.)

z %>% 
  filter(gene == "asip")

z %>% 
  pivot_longer(corin:ednrb, names_to = "gene", values_to = "dosage") %>% 
  pivot_longer(snow.cover:bio_3, names_to = "env_var", values_to = "value") %>% 
  ggplot(aes(x = value, y = dosage, color = gene)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_grid(cols = vars(env_var), scales = "free_x")


z %>% 
  pivot_longer(corin:ednrb, names_to = "gene", values_to = "dosage") %>% 
  pivot_longer(snow.cover:bio_3, names_to = "env_var", values_to = "value") %>% 
  filter(gene == "asip") %>% 
  filter(env_var == "snow.cover")


cor.test(z$snow.cover, z$corin)
cor.test(z$bio_2, z$corin)
cor.test(z$bio_3, z$corin)


cor.test(z$snow.cover, z$asip)
cor.test(z$bio_2, z$asip)
cor.test(z$bio_3, z$asip)


cor.test(z$snow.cover, z$ednrb)
cor.test(z$bio_2, z$ednrb)
cor.test(z$bio_3, z$ednrb)

tidy(lm(corin ~ snow.cover + bio_2 + bio_3, data = z))
  

env_sites <- z %>% 
  dplyr::select(latitude, longitude, snow.cover:bio_3) %>% 
  distinct() %>% 
  mutate(site = 1:nrow(.))

env_pca <- env_sites %>% 
  dplyr::select(-latitude, -longitude) %>% 
  column_to_rownames("site") %>%
  prcomp(., scale = T)
  
env_pca$x

env_sites_w_pcs <- bind_cols(env_sites, as.data.frame(env_pca$x))
  
combo_2 <- wtjr_locs %>% 
  left_join(env_sites_w_pcs)



combo_2 %>% 
  dplyr::select(corin:ednrb, PC1:PC3) %>% 
  ggpairs(.)


tidy(lm(corin ~ snow.cover + bio_2 + bio_3, data = z))
tidy(lm(cbind(corin, asip, ednrb) ~ snow.cover + bio_2 + bio_3, data = z)) %>% 
  mutate(sig = p.value < 0.05)




fviz_pca_biplot(env_pca, repel = T) + theme_cowplot()


# Preliinaryilty, seems promising.
# Genotupe dosages at all 3 genes are relativley strongly correlated with snow cover
# interestingly, snow cover is not that strongly correlated with latitutde across our sites, though it is a
# little negatively correlated with longitude.

# Patterns at bioclim 2 and bioclim 3 don't make as much sense. Associations can be strong, but 
# the directionality doesn't always make a lot of sense given the directionality of the beta coefficients
# in the predictive phenotype models. 
# But this could be because bioclim 2 and 3 are approaching colinearity across our environmental range
# and both are approaching colinearity with latitude. 

  


# Figure out if PCAngsd dosages are for major or minor allele. 

gl_sample_order <- paste0(unlist(lapply(seqed_samples$sample, rep, times = 3)), c("_major", "_het", "_minor"))
asip_gls <- read_delim("results/pcangsd/all_wtjr_samples/asip_gl_line.txt", delim = "\t", 
                       col_names = c("position", "allele1", "allele2", gl_sample_order))

z <- asip_gls %>% 
  pivot_longer(4:ncol(.)) 

%>% 
  filter(str_detect(name, "DMNS-4313"))
