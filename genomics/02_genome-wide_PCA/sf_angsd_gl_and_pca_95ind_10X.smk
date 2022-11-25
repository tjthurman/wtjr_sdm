# Snakemake pipeline to get a PCA of the WJTR BAMs with ANGSD
# After rule "all", pipeline order goes from top to bottom

###############
##   SETUP   ##
############### 
import os
import re
# import csv


###############
##  GLOBALS  ##
###############

# Reference genome
REF="/mnt/beegfs/tt164677e/genomes/lepus_townsendii/DMNS18807_06042020_pseudohap2.1_10k.fasta"

# Get scaffold list:
with open(REF + ".fai", "r") as f:
    lines =  f.read().splitlines()
    scaffolds = [line.split("\t")[0] for line in lines]

# For testing on a small number of scafs
# scaffolds = scaffolds[30:36]

####################
## PIPELINE START ##
####################
localrules: all


rule all:
    input:
        "results/angsd_95ind/beagle_GL_scaffold/all_wtjr_samples/genome_wide_GL_all_samples.beagle.gz",
        "results/pcangsd_95ind/all_wtjr_samples/pca_genome_wide_GL_all_samples.cov"
       


# calculate genotype likelihoods separately for each scaffold
# Simultanesouly does filering,
# calculating of genotype 
rule calc_beagle_GLs_by_scaff:
    input:
        bamlist= "data/bamlists/all_samples_bamlist.txt",
        ref=REF
    output:
        "results/angsd_95ind/beagle_GL_scaffold/all_wtjr_samples/by_scaffold/{scaffold}_GL.beagle.gz"
    params:
        basename=lambda wildcards: "results/angsd_95ind/beagle_GL_scaffold/all_wtjr_samples/by_scaffold/" + wildcards.scaffold + "_GL",
        region=lambda wildcards: wildcards.scaffold + ":"
    priority: 10
    log:
        "logs/angsd_GL_by_scaff/{scaffold}.txt"
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
            -minInd 95 -setMinDepth 95 -setMaxDepth 1430 \
            -nThreads  {resources.cpus}
        """

rule gather_per_scaff_GLS:
    input:
        expand("results/angsd_95ind/beagle_GL_scaffold/all_wtjr_samples/by_scaffold/{scaffold}_GL.beagle.gz", scaffold = scaffolds)
    output:
        glf="results/angsd_95ind/beagle_GL_scaffold/all_wtjr_samples/genome_wide_GL_all_samples.beagle.gz",
        sites="results/angsd_95ind/beagle_GL_scaffold/all_wtjr_samples/genome_wide_GL_all_samples.beagle.sites"
    shell:
        """
        cat {input} > {output.glf}
        zcat {output.glf} | awk -F '\t' '{{print $1}}' > {output.sites}
        """


# Using the wtjr_analyze_bams conda environment for this rule
rule pcangsd:
    input:
        "results/angsd_95ind/beagle_GL_scaffold/all_wtjr_samples/genome_wide_GL_all_samples.beagle.gz"
    output:
        "results/pcangsd_95ind/all_wtjr_samples/pca_genome_wide_GL_all_samples.cov"
    params:
        basename="results/pcangsd_95ind/all_wtjr_samples/pca_genome_wide_GL_all_samples"
    resources:
        cpus=36,
        partition="good_lab_reincarnation",
        mem_mb=100000
    shell:
        """
        python bin/pcangsd/pcangsd.py -beagle {input} -out {params.basename} -threads {resources.cpus} -pcadapt -selection -sites_save -dosage_save
        """