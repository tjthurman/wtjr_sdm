import numpy as np
import pandas as pd

dosages = np.load("results/pcangsd/all_wtjr_samples/pca_genome_wide_GL_all_samples.dosage.npy", mmap_mode='r')

# So, this memory-maps an array but doesn't read it on to disk, so I don't have to worry about reading 
# a file that is 14Gb when compressed

dosages.shape
# Give me the dimensions
# it is 25779497 rows (25 million). Thats the number of SNPs it calculated dosages for
# it is 143 columns. Those are the individuals. 
# PCANGSD doesn't explicitly list the individual order.
# But, I assume it is the same as the bamlist. In any case, that is something to worry about later.
# First, need to figure out the positions of each of the sites

# So, we need to look at the sites file output by PCANGSD
# It outputs a either a 1 or a 0, for whether a site was included or not. 
pcangsd_sites = np.loadtxt("results/pcangsd/all_wtjr_samples/pca_genome_wide_GL_all_samples.sites")
pcangsd_sites.shape
# so, this file is 41814111 sites (41 million). 
# The sum of this should be the same as the number of rows in the dosages
np.sum(pcangsd_sites)
# And it is. So, that's good. 

# I can also read it in as text
with open("results/pcangsd/all_wtjr_samples/pca_genome_wide_GL_all_samples.sites") as f:
    pcangsd_sites_text = f.read().splitlines()
# And it is the same size


# But this file just lists if a site was included or not. 
# it doesn't give the contig and position
# That is in the first column of the beagle.gz file
# And also saved as the .beagle.sites file form the angsd output

# I can load the .beagle.sites file
with open("results/angsd/beagle_GL_scaffold/all_wtjr_samples/genome_wide_GL_all_samples.beagle.sites") as f:
    all_positions = f.read().splitlines()
# It has 41814112 sites listed in it
# This is one more than the list of sites output by PCANGSD. 
# In spot-checking the the actual text files against each other, there is 
# definitely some sort of off-by-one error SOMEWHERE, as some sites named "marker"
# are included, even though that is the head of the filename. 

# Let's try using zcat and awk to pull the file markers
# directly from the beagle.gz file that went into PCANGSD.
# This might take a while...
zcat genome_wide_GL_all_samples.beagle.gz | awk '{print $1}' > sites_into_pcangsd.txt

# This didn't work
# Still one extra line
# Let's try to find the extra line.
# The input beagle.gz contains rows of "marker", which are the headers of the individual beagle.gz files that were concatenated
# together to make the genome-wide file
# Those sites that are "marker" should be excluded by PCANGSD

# Heres a horrible awk serach that finds the indices of all the marker sites in the output from angsd
# and then returns the corresponding 0/1 value from PCANGSD. It shoudl return 1: finding the first spot where
# they don't should help me diagnose where the issue occurs. 
awk 'NR==FNR { if ($1 == "marker") {out[NR]= 1; next}}  { if (out[FNR]==1) {print (FNR, $0);}  }'  results/angsd/beagle_GL_scaffold/all_wtjr_samples/genome_wide_GL_all_samples.beagle.sites results/pcangsd/all_wtjr_samples/pca_genome_wide_GL_all_samples.sites > awk_search_results.txt

# Looks like it is basically after the first scaffold with any data in it (84). 
# Seems like the issue might be with the very first line????
# Can remove it:
sed '1d' sites_into_pcangsd.txt > sites_into_pcangsd_fixed.txt

# Let's remove it and see if all the markers have 0 agian. 
awk 'NR==FNR { if ($1 == "marker") {out[NR]= 1; next}}  { if (out[FNR]==1) {print (FNR, $0);}  }'  results/angsd/beagle_GL_scaffold/all_wtjr_samples/sites_into_pcangsd_fixed.txt results/pcangsd/all_wtjr_samples/pca_genome_wide_GL_all_samples.sites > awk_search_results_2.txt
# They are all 0! Seems to have fixed it

# So, I can re-load the beagle positions
with open("results/angsd/beagle_GL_scaffold/all_wtjr_samples/sites_into_pcangsd_fixed.txt") as f:
    all_positions = f.read().splitlines()
# And that gets the right length


# Now, to finally do the subsetting:
# First, find where the pcangsd output is equal to 1
# do some dumb selection to get it to not be a tuple
keepers = np.where(pcangsd_sites == 1)[0]

# and a list comprehension to turn it into a list
kept_positions = [all_positions[i] for i in keepers.tolist()]

# Then, find the index where the positions cooridinates
# match those of the top associated SNPs for each gene
corin_sites = ["342_47124004","342_46991691","342_46991393","342_46989842","342_46983691","342_46972767","342_46970080","342_46966737","342_46958289","342_46957711","342_46879687","342_46871251","342_46870994"]

# Not found in list of sites: 245_24228244
asip_sites = ["245_24236852","245_24229551","245_24226798","245_24225226","245_24222426","245_24221492","245_24220190","245_24211740","245_24204540","245_24200684","245_24197147"]

# Not found in sites: 311_3456181
EDNRB_sites = ["311_3467105","311_3456263","311_3454568","311_3419309","311_3419220","311_3413288","311_3402184","311_3399522"]

# Not found in site: 380_35303804, 380_35322019
unk_sites = ["380_35303819", "380_35321971", "380_35322089", "380_35311674", "380_35320825", "380_35253666", "380_35313330", "380_35306181"]


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
corin_dosage_df.to_csv("results/pcangsd/all_wtjr_samples/corin_dosages.csv", header = False)
asip_dosage_df.to_csv("results/pcangsd/all_wtjr_samples/asip_dosages.csv", header = False)
ednrb_dosage_df.to_csv("results/pcangsd/all_wtjr_samples/ednrb_dosages.csv", header = False)
unk_dosage_df.to_csv("results/pcangsd/all_wtjr_samples/scaff380_dosages.csv", header = False)




