
library(tidyverse)

# Get samples sequenced ---------------------------------------------------
# Import from the bamlist
seqed_samples <- read_delim(file = "data/all_wtjr_sample_bamlist.txt", delim = "\t", col_names = "bam_file") %>% 
  mutate(sample = basename(bam_file)) %>% 
  mutate(sample = str_remove(sample, ".realigned.sorted.bam")) %>% 
  mutate(sample = str_remove(sample, "_realigned_sorted.bam"))

samples <- seqed_samples %>% 
  pull(sample)  



# Import beagle GLs -------------------------------------------------------
# Also gives us allele identity
# A function to convert the numric allel identites from angsd to characters
angsd_allele_ids <- function(numeric) {
  switch (numeric,
    "0" = {out <- "A"},
    "1" = {out <- "C"},
    "2" = {out <- "G"},
    "3" = {out <- "T"},
    stop("Not a valid entry")
  )
  out
}
angsd_allele_ids <- Vectorize(angsd_allele_ids)
angsd_allele_ids("0")

  
# Import each
asip_GLs <- read_delim("results/angsd/beagle_GL_scaffold/all_wtjr_samples/asip_site_GLs.beagle", delim = "\t", 
           col_names = c("location", 
                         "major_allele",
                         "minor_allele",
                         paste(rep(samples, each = 3), c("MM", "Mm", "mm"), sep = "_"))) %>% 
  separate(location, into = c("scaffold", "position"), remove = F, convert = T) %>% 
  mutate(major_allele = angsd_allele_ids(as.character(major_allele))) %>% 
  mutate(minor_allele = angsd_allele_ids(as.character(minor_allele)))


corin_GLs <- read_delim("results/angsd/beagle_GL_scaffold/all_wtjr_samples/corin_site_GLs.beagle", delim = "\t", 
                       col_names = c("location", 
                                     "major_allele",
                                     "minor_allele",
                                     paste(rep(samples, each = 3), c("MM", "Mm", "mm"), sep = "_"))) %>% 
  separate(location, into = c("scaffold", "position"), remove = F, convert = T) %>% 
  mutate(major_allele = angsd_allele_ids(as.character(major_allele))) %>% 
  mutate(minor_allele = angsd_allele_ids(as.character(minor_allele))) 

ednrb_GLs <- read_delim("results/angsd/beagle_GL_scaffold/all_wtjr_samples/ednrb_site_GLs.beagle", delim = "\t", 
                       col_names = c("location", 
                                     "major_allele",
                                     "minor_allele",
                                     paste(rep(samples, each = 3), c("MM", "Mm", "mm"), sep = "_"))) %>% 
  separate(location, into = c("scaffold", "position"), remove = F, convert = T) %>% 
  mutate(major_allele = angsd_allele_ids(as.character(major_allele))) %>% 
  mutate(minor_allele = angsd_allele_ids(as.character(minor_allele))) 

scaff380_GLs <- read_delim("results/angsd/beagle_GL_scaffold/all_wtjr_samples/scaff380_site_GLs.beagle", delim = "\t", 
                           col_names = c("location", 
                                         "major_allele",
                                         "minor_allele",
                                         paste(rep(samples, each = 3), c("MM", "Mm", "mm"), sep = "_"))) %>% 
  separate(location, into = c("scaffold", "position"), remove = F, convert = T) %>% 
  mutate(major_allele = angsd_allele_ids(as.character(major_allele))) %>% 
  mutate(minor_allele = angsd_allele_ids(as.character(minor_allele))) 



# Figure out which site is most imformative -------------------------------
# Might be interesting to try to think of this in terms of bits,
# just out of interest.
# Then, could take the site with the most bits. 

# Perhaps too much to worry about for now. but if I can figure out an easy way, maybe. 
asip_GLs %>% 
  dplyr::select(-(scaffold:minor_allele)) %>% 
  pivot_longer(-location, names_to =  "sample") %>% 
  mutate(sample = str_remove(sample, "_[[Mm]]+")) %>% 
  mutate(dist = (value - (1/3))^2) %>% 
  group_by(location, sample) %>% 
  summarize(qc_score = sqrt(sum(dist))/0.8164966) %>% 
  ungroup() %>% 
  write_csv("results/angsd/asip_GL_qc_scores.csv")


corin_GLs %>% 
  dplyr::select(-(scaffold:minor_allele)) %>% 
  pivot_longer(-location, names_to =  "sample") %>% 
  mutate(sample = str_remove(sample, "_[[Mm]]+")) %>% 
  mutate(dist = (value - (1/3))^2) %>% 
  group_by(location, sample) %>% 
  summarize(qc_score = sqrt(sum(dist))/0.8164966) %>% 
  ungroup() %>% 
  write_csv("results/angsd/corin_GL_qc_scores.csv")

ednrb_GLs %>% 
  dplyr::select(-(scaffold:minor_allele)) %>% 
  pivot_longer(-location, names_to =  "sample") %>% 
  mutate(sample = str_remove(sample, "_[[Mm]]+")) %>% 
  mutate(dist = (value - (1/3))^2) %>% 
  group_by(location, sample) %>% 
  summarize(qc_score = sqrt(sum(dist))/0.8164966) %>% 
  ungroup() %>% 
  write_csv("results/angsd/ednrb_GL_qc_scores.csv")

scaff380_GLs %>% 
  dplyr::select(-(scaffold:minor_allele)) %>% 
  pivot_longer(-location, names_to =  "sample") %>% 
  mutate(sample = str_remove(sample, "_[[Mm]]+")) %>% 
  mutate(dist = (value - (1/3))^2) %>% 
  group_by(location, sample) %>% 
  summarize(qc_score = sqrt(sum(dist))/0.8164966) %>% 
  ungroup() %>% 
  write_csv("results/angsd/scaff380_GL_qc_scores.csv")
  



# Import Massarray info on allele identity --------------------------------
# Massarry is based on orycun coords,
# which are the reverse complement of WTJR
# So, need to complement them here


complement <- function(dnachar) {
  switch(dnachar,
         "A" = out <- "T",
         "C" = out <- "G",
         "G" = out <- "C",
         "T" = out <- "A",
         stop("Not valid base")
         )
  out
}
complement <- Vectorize(complement)

allele_ids <- read_csv("data/white_brown_allele_ids.csv") %>% 
  filter(!is.na(gene)) %>% 
  mutate(white_allele = complement(white_allele)) %>% 
  mutate(brown_allele = complement(brown_allele))


# Import info on translating between coord systems ------------------------
translate_coords <- read_csv("data/color_allele_coord_translation.csv", col_types = "cc") %>% 
  separate(orycun_location, into = c("scaffold_orycun", "position_orycun"), remove = T, convert = T) %>% 
  separate(wtjr_location, into = c("scaffold", "position"), remove = T, convert = T)


# Get "color" of minor allele for each site -------------------------------
# Sort of an ungodly way to do it, but works for now
asip_colors <- asip_GLs %>% 
  select(location, scaffold, position, major_allele, minor_allele) %>% 
  left_join(translate_coords) %>% 
  left_join(allele_ids) %>% 
  mutate(brown = minor_allele == brown_allele) %>% 
  mutate(white = minor_allele == white_allele) %>% 
  mutate(other = minor_allele != brown_allele & minor_allele != white_allele) %>% 
  dplyr::select(location, scaffold, position, major_allele, brown, white, other) %>% 
  pivot_longer(brown:other, names_to = "minor_color") %>% 
  filter(value) %>% 
  select(-value) %>% 
  left_join(translate_coords) %>% 
  left_join(allele_ids) %>%
  mutate(brown = major_allele == brown_allele) %>% 
  mutate(white = major_allele == white_allele) %>% 
  mutate(other = major_allele != brown_allele & major_allele != white_allele) %>% 
  dplyr::select(location, scaffold, position, minor_color, brown:other) %>% 
  pivot_longer(brown:other, names_to = "major_color") %>% 
  filter(value) %>% 
  select(-value) 

corin_colors <- corin_GLs %>% 
  select(location, scaffold, position, major_allele, minor_allele) %>% 
  left_join(translate_coords) %>% 
  left_join(allele_ids) %>% 
  mutate(brown = minor_allele == brown_allele) %>% 
  mutate(white = minor_allele == white_allele) %>% 
  mutate(other = minor_allele != brown_allele & minor_allele != white_allele) %>% 
  dplyr::select(location, scaffold, position, major_allele, brown, white, other) %>% 
  pivot_longer(brown:other, names_to = "minor_color") %>% 
  filter(value) %>% 
  select(-value) %>% 
  left_join(translate_coords) %>% 
  left_join(allele_ids) %>%
  mutate(brown = major_allele == brown_allele) %>% 
  mutate(white = major_allele == white_allele) %>% 
  mutate(other = major_allele != brown_allele & major_allele != white_allele) %>% 
  dplyr::select(location, scaffold, position, minor_color, brown:other) %>% 
  pivot_longer(brown:other, names_to = "major_color") %>% 
  filter(value) %>% 
  select(-value) 


ednrb_colors <- ednrb_GLs %>% 
  select(location, scaffold, position, major_allele, minor_allele) %>% 
  left_join(translate_coords) %>% 
  left_join(allele_ids) %>% 
  mutate(brown = minor_allele == brown_allele) %>% 
  mutate(white = minor_allele == white_allele) %>% 
  mutate(other = minor_allele != brown_allele & minor_allele != white_allele) %>% 
  dplyr::select(location, scaffold, position, major_allele, brown, white, other) %>% 
  pivot_longer(brown:other, names_to = "minor_color") %>% 
  filter(value) %>% 
  select(-value) %>% 
  left_join(translate_coords) %>% 
  left_join(allele_ids) %>%
  mutate(brown = major_allele == brown_allele) %>% 
  mutate(white = major_allele == white_allele) %>% 
  mutate(other = major_allele != brown_allele & major_allele != white_allele) %>% 
  dplyr::select(location, scaffold, position, minor_color, brown:other) %>% 
  pivot_longer(brown:other, names_to = "major_color") %>% 
  filter(value) %>% 
  select(-value) 


# With scaffold 380, don't have any preset colors. Need to figure this out. 
# First, will read in the QC scores and subset to only sites/individuals with good enough QC scores to try to 
# figure out which allele is brown and which is white. 
scaff380_qc <- read_csv("results/angsd/scaff380_GL_qc_scores.csv")
# Then, pull in the color info for each sample from table S3
sample_colors <- read_csv("data/sample_colors.csv") %>% 
  mutate(sample = str_replace(sample, "_", "-"))
x <- scaff380_GLs %>% 
  pivot_longer(AMNH_123864_MM:last_col()) %>% 
  mutate(sample = str_remove(name, "_[[Mm]]+")) %>% 
  mutate(genotype = str_extract(name, "_[[Mm]]+")) %>% 
  left_join(scaff380_qc) %>% 
  filter(qc_score > 0.6) %>% 
  left_join(sample_colors) %>% 
  filter(!is.na(color)) %>% 
  group_by(location, scaffold, position, major_allele, minor_allele, genotype, color) %>% 
  summarize(mean_gl = mean(value, na.rm = T))


scaff380_dosages <- read_csv("results/pcangsd/all_wtjr_samples/scaff380_dosages.csv", col_names = c("location", samples)) %>% 
  separate(location, into = c("scaffold", "position"), remove = F, convert = T) 


z <- scaff380_dosages %>% 
  pivot_longer(AMNH_123864:last_col(), names_to = "sample", values_to = "dosage") %>% 
  left_join(sample_colors) %>% 
  filter(!is.na(color)) %>% 
  group_by(location, scaffold, position, color) %>% 
  summarize(mean_dosage = mean(dosage))




# Import Genotype dosages -------------------------------------------------
asip_dosages <- read_csv("results/pcangsd/all_wtjr_samples/asip_dosages.csv", col_names = c("location", samples)) %>% 
  separate(location, into = c("scaffold", "position"), remove = F, convert = T) %>% 
  left_join(asip_colors) %>% 
  relocate(minor_color:major_color, .after = position)

corin_dosages <- read_csv("results/pcangsd/all_wtjr_samples/corin_dosages.csv", col_names = c("location", samples)) %>% 
  separate(location, into = c("scaffold", "position"), remove = F, convert = T) %>% 
  left_join(corin_colors) %>% 
  relocate(minor_color:major_color, .after = position)

ednrb_dosages <- read_csv("results/pcangsd/all_wtjr_samples/ednrb_dosages.csv", col_names = c("location", samples)) %>% 
  separate(location, into = c("scaffold", "position"), remove = F, convert = T) %>% 
  left_join(ednrb_colors) %>% 
  relocate(minor_color:major_color, .after = position)


# Get all dosages in terms of white
dosages_white <- function(dosage_row) {
  if (dosage_row["minor_color"] == "white") { # if minor is already white
    return(dosage_row) # don't do anything
  } else if (dosage_row["major_color"] == "white") { # if major is white, switch dosage 
    dosage_white <- dosage_row
    dosage_white[6:length(dosage_row)] <- 2 - as.numeric(dosage_row[6:length(dosage_row)])
    return(dosage_white)
  } else { # minor was other and major was brown
    # Can't use a site like that
    dosage_out <- dosage_row
    dosage_out[6:length(dosage_row)] <- NA
    return(dosage_out)
  }
}


asip_white_dosage <- as_tibble(t(apply(X = asip_dosages, MARGIN = 1, FUN = dosages_white)))
corin_white_dosage <- as_tibble(t(apply(X= corin_dosages, MARGIN = 1, FUN = dosages_white)))
ednrb_white_dosage <- as_tibble(t(apply(X= ednrb_dosages, MARGIN = 1, FUN = dosages_white)))

scaff380_white_dosage <- scaff380_dosages %>% 
  pivot_longer(AMNH_123864:last_col()) %>% 
  mutate(value = 2 - value) %>% 
  pivot_wider()

# Save them all
asip_white_dosage %>% 
  select(-minor_color, -major_color) %>% 
  write_csv("results/pcangsd/asip_dosages_white.csv", col_names = T)

corin_white_dosage %>% 
  select(-minor_color, -major_color) %>% 
  write_csv("results/pcangsd/corin_dosages_white.csv", col_names = T)

ednrb_white_dosage %>% 
  select(-minor_color, -major_color) %>% 
  write_csv("results/pcangsd/ednrb_dosages_white.csv", col_names = T)

scaff380_white_dosage %>% 
  write_csv("results/pcangsd/scaff380_dosages_white.csv", col_names = T)


# After working by hand, seems like scaff380's minor alleles are all for the brown allele, so will just do 2- dosage to save

