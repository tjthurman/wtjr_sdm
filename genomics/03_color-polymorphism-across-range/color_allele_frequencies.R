library(tidyverse)
library(ggrepel)


# Import data on allele identify ------------------------------------------
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

allele_ids <- read_csv("genomics/03_color-polymorphism-across-range/white_brown_allele_ids.csv") %>% 
  filter(!is.na(gene)) %>% 
  mutate(white_allele = complement(white_allele)) %>% 
  mutate(brown_allele = complement(brown_allele))


# Import info on translating between coord systems ------------------------
translate_coords <- read_csv("genomics/03_color-polymorphism-across-range/color_allele_coord_translation.csv", col_types = "cc") %>% 
  separate(orycun_location, into = c("scaffold_orycun", "position_orycun"), remove = T, convert = T) %>% 
  separate(wtjr_location, into = c("scaffold", "position"), remove = T, convert = T)


# Load in allele frequency data -------------------------------------------

# ASIP --------------------------------------------------------------------
top_asip <- allele_ids %>% 
  left_join(translate_coords) %>% 
  filter(position == 24236852)


asip_nonCO <- read_delim("genomics/results/angsd_by_region/not_colorado/asip_GL.mafs.gz", delim = "\t") %>% 
  filter(position == 24236852) %>% 
  left_join(top_asip) %>% 
  mutate(brown_freq = unknownEM,
         white_freq = 1-brown_freq,
         gene = "asip",
         pop = "outside CO") %>% 
  select(gene, pop, white_freq, brown_freq)

asip_CO <- read_delim("genomics/results/angsd_by_region/colorado/asip_GL.mafs.gz", delim = "\t") %>% 
  filter(position == 24236852) %>% 
  left_join(top_asip) %>% 
  mutate(brown_freq = unknownEM,
         white_freq = 1-brown_freq,
         gene = "asip",
         pop = "CO") %>% 
  select(gene, pop, white_freq, brown_freq)



# Corin --------------------------------------------------------------------
top_corin <- allele_ids %>% 
  left_join(translate_coords) %>% 
  filter(position == 46966737)


corin_nonCO <- read_delim("genomics/results/angsd_by_region/not_colorado/corin_GL.mafs.gz", delim = "\t") %>% 
  filter(position == 46966737) %>% 
  left_join(top_corin) %>% 
  mutate(white_freq = unknownEM,
         brown_freq = 1-white_freq,
         gene = "corin",
         pop = "outside CO") %>% 
  select(gene, pop, white_freq, brown_freq)

corin_CO <- read_delim("genomics/results/angsd_by_region/colorado/corin_GL.mafs.gz", delim = "\t") %>% 
  filter(position == 46966737) %>% 
  left_join(top_corin) %>% 
  mutate(white_freq = unknownEM,
         brown_freq = 1-white_freq,
         gene = "corin",
         pop = "CO") %>% 
  select(gene, pop, white_freq, brown_freq)




# EDNRB -------------------------------------------------------------------

top_ednrb <- allele_ids %>% 
  left_join(translate_coords) %>% 
  filter(position == 3413288)


ednrb_nonCO <- read_delim("genomics/results/angsd_by_region/not_colorado/ednrb_GL.mafs.gz", delim = "\t") %>% 
  filter(position == 3413288) %>% 
  left_join(top_ednrb) %>% 
  mutate(brown_freq = unknownEM,
         white_freq = 1-brown_freq,
         gene = "ednrb",
         pop = "outside CO") %>% 
  select(gene, pop, white_freq, brown_freq)

ednrb_CO <- read_delim("genomics/results/angsd_by_region/colorado/ednrb_GL.mafs.gz", delim = "\t") %>% 
  filter(position == 3413288) %>% 
  left_join(top_ednrb) %>% 
  mutate(white_freq = unknownEM,
         brown_freq = 1-white_freq,
         gene = "ednrb",
         pop = "CO") %>% 
  select(gene, pop, white_freq, brown_freq)



# Top scaff 380 -----------------------------------------------------------
# Need to figure out white allele and brown allele

top_scaff380 <- allele_ids %>% 
  left_join(translate_coords) %>% 
  filter(position == 35303819)

# In 380 at the top site, the minor allele, G, is associated with B individuals
# in the few individuals we have info for.
# This goes along with the dosages in general: pcangsd gives the does of the minor allele,
# and for the overall sample the minor allele was G.
# In these data, both CO and nonCO the minor allele is G, so its the brown allele. 

scaff380_nonCO <- read_delim("genomics/results/angsd_by_region/not_colorado/scaff380_GL.mafs.gz", delim = "\t") %>% 
  filter(position == 35303819) %>% 
  mutate(brown_freq = unknownEM,
         white_freq = 1-brown_freq,
         gene = "scaff380",
         pop = "outside CO") %>% 
  select(gene, pop, brown_freq, white_freq)

scaff380_CO <- read_delim("genomics/results/angsd_by_region/colorado/scaff380_GL.mafs.gz", delim = "\t") %>% 
  filter(position == 35303819) %>% 
  mutate(brown_freq = unknownEM,
         white_freq = 1-brown_freq,
         gene = "scaff380",
         pop = "CO") %>% 
  select(gene, pop, brown_freq, white_freq)



for_pie <- rbind(asip_nonCO, asip_CO, 
      corin_nonCO, corin_CO, 
      ednrb_nonCO, ednrb_CO,
      scaff380_nonCO, scaff380_CO) %>%
  pivot_longer(white_freq:brown_freq, names_to = "color") %>% 
  mutate(color = str_remove(color, "_freq")) %>% 
  mutate(region = ifelse(pop == "CO", "Colorado", "outside Colorado")) %>% 
  mutate(gene = ifelse(gene == "scaff380", gene, str_to_upper(gene))) 

label_pos <- for_pie %>%
  group_by(gene, pop) %>% 
  mutate(pos = value/2 + lag(value, 1),
         pos = if_else(is.na(pos), value/2, pos))


for_pie %>% 
  ggplot(aes(x = " ", y = value)) +
  geom_bar(aes( group = color, fill = color), width = 1, stat = "identity", color = "black") +
  coord_polar("y", start=0) + 
  facet_grid(gene ~ region) + 
  geom_text(data = label_pos,
             aes(y = pos,
                 label = paste0(as.integer(value*100), "%"),
                 color = color), 
             show.legend = FALSE) +
  scale_fill_manual(values = c("brown" = "#69431A",
                               "white" = "floralwhite")) +
  scale_color_manual(values = c("brown" = "white",
                                "white" = "black")) +
  theme_void() + 
  theme(legend.position = "none",
        strip.text = element_text(size = 12))
ggsave("genomics/results/figures/supplemental/freq_piecharts.pdf", width = 4, height = 11, units = "in")
