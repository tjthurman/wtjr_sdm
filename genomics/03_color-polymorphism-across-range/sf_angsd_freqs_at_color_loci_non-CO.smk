# Snakemake pipeline to calculate allele frequencies just in our color regions
# for samples outside CO
# run from angsd conda environment. 

###############
##  GLOBALS  ##
###############

# Reference genome
REF="/mnt/beegfs/tt164677e/genomes/lepus_townsendii/DMNS18807_06042020_pseudohap2.1_10k.fasta"

# Set up dictionary of regions in which to calculate allele freqs
# Just the color-associated regions

genes = {
    "asip": "245:24197146-24236853",
    "corin": "342:46879686-47124005",
    "ednrb": "311:3399521-3467106",
    "scaff380": "380:35303803-35322020"
}

####################
## PIPELINE START ##
####################
localrules: all


rule all:
    input:
        expand("results/angsd_by_region/not_colorado/{gene}_GL.beagle.gz", gene = genes.keys())
       
# calculate genotype likelihoods separately for each scaffold
# Simultanesouly does filering,
# calculating of genotype 
rule calc_beagle_GLs_by_scaff:
    input:
        bamlist= "data/bamlists/non-colorado_samples_bamlist.txt",
        ref=REF
    output:
        "results/angsd_by_region/not_colorado/{gene}_GL.beagle.gz"
    params:
        basename=lambda wildcards: "results/angsd_by_region/not_colorado/" + wildcards.gene + "_GL",
        region=lambda wildcards: genes.get(wildcards.gene)
    priority: 10
    log:
        "logs/angsd_by_region/not_colorado/{gene}.txt"
    resources:
        cpus=18
    shell:
        """
        angsd -bam {input.bamlist} \
            -ref {input.ref} \
            -anc {input.ref} \
            -out {params.basename} \
            -r {params.region} \
            -gl 1 -doGlf 2 -doCounts 1 -dosaf 1 \
            -doMajorMinor 1 -doMaf 2 \
            -SNP_pval 1e-6 \
            -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
            -trim 0 -C 50 -baq 1 \
            -nThreads  {resources.cpus}
        """
