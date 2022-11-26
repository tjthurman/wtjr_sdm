import numpy as np
import pandas as pd




# This file is too large to store in Git, and is output by PCAngsd
dosages = np.load("results/pcangsd_95ind/all_wtjr_samples/pca_genome_wide_GL_all_samples.dosage.npy", mmap_mode='r')
dosages.shape # Number of SNPs in the PCA analysis. 

# Lists which sites were and weren't included in the PCA analysis
# This file is also output by PCAngsd. 0 indicates site wasn't include, 1 indicates it was.
pcangsd_sites = np.loadtxt("results/pcangsd_95ind/all_wtjr_samples/pca_genome_wide_GL_all_samples.sites")
pcangsd_sites.shape


# Then, load in the position IDs of all sites that went into PCANGSD,
# output from the get_pcangsd_sites.sh script. 
with open("results/angsd_95ind/beagle_GL_scaffold/all_wtjr_samples/sites_into_pcangsd_fixed.txt") as f:
    all_positions = f.read().splitlines()


# On to the subsetting:
# First, find where the pcangsd output is equal to 1
keepers = np.where(pcangsd_sites == 1)[0]

# and a list comprehension to turn it into a list
kept_positions = [all_positions[i] for i in keepers.tolist()]



# Then, find the index where the positions cooridinates
# match those of the top associated SNPs for each gene
corin_sites = ["342_47124004", "342_46991393", "342_46989842", "342_46983691", "342_46972767", "342_46966737", "342_46958289", "342_46957711", "342_46879687", "342_46871251","342_46870994"]
asip_sites = [ "245_24236852","245_24229551","245_24225226","245_24221492","245_24211740","245_24204540","245_24197147"]
EDNRB_sites = [ "311_3467105","311_3456263","311_3454568","311_3419309","311_3419220","311_3413288","311_3402184","311_3399522"]
unk_sites = [ "380_35303819", "380_35322089", "380_35311674", "380_35320825", "380_35253666", "380_35313330", "380_35306181"]


corin_indices = [kept_positions.index(x) for x in corin_sites]
ASIP_indices = [kept_positions.index(x) for x in asip_sites]
EDNRB_indices = [kept_positions.index(x) for x in EDNRB_sites]
unk_indices = [kept_positions.index(x) for x in unk_sites]

# Then, pull out the dosages for these
corin_dosages = dosages[corin_indices, :]
asip_dosages = dosages[ASIP_indices, :]
ednrb_dosages = dosages[EDNRB_indices, :]
unk_dosages = dosages[unk_indices, :]


# Convert to data frames
corin_dosage_df = pd.DataFrame(corin_dosages, index = corin_sites)
asip_dosage_df = pd.DataFrame(asip_dosages, index = asip_sites)
ednrb_dosage_df = pd.DataFrame(ednrb_dosages, index = EDNRB_sites)
unk_dosage_df = pd.DataFrame(unk_dosages, index = unk_sites)

# And save them as tsv 
corin_dosage_df.to_csv("genomics/results/pcangsd/95ind_10X_filter/corin_dosages.csv", header = False)
asip_dosage_df.to_csv("genomics/results/pcangsd/95ind_10X_filter/asip_dosages.csv", header = False)
ednrb_dosage_df.to_csv("genomics/results/pcangsd/95ind_10X_filter/ednrb_dosages.csv", header = False)
unk_dosage_df.to_csv("genomics/results/pcangsd/95ind_10X_filter/scaff380_dosages.csv", header = False)

