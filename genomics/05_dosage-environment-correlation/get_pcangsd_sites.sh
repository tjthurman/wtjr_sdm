

# Change to the approtpariate directory that contains the
# genome_wide_GL_all_samples.beagle.gz file output by the snakemake
# pipline in the `02_genome-wide_PCA` folder



# From the genotype likelihood file of all sites that went into PCANGSD, extract the 
# IDs of the sites
zcat genome_wide_GL_all_samples.beagle.gz | awk '{print $1}' > sites_into_pcangsd.txt

# And fix an issue with the first line so that line IDs match up properly later. 
sed '1d' sites_into_pcangsd.txt > sites_into_pcangsd_fixed.txt