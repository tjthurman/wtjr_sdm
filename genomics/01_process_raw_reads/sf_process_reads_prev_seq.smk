# Initial Snakemake pipeline to process reads data in prep for GEAs
# After rule "all", pipeline order goes from top to bottom

# Downloaded old GATK 3.6 from
# https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

###############
##   SETUP   ##
############### 
localrules: all, help, dag

# Common variables:
# Reference genome
REF="/mnt/beegfs/tt164677e/genomes/lepus_townsendii/DMNS18807_06042020_pseudohap2.1_10k.fasta"

# Path to genome index file for BWA
index_path= REF + ".amb"

# Samples inputs and lists:

# Python list of samples with their lane and run info
# This will serve as the list of targets for pre-merging steps
with open("processed_data/sample_with_lane_and_run.txt") as f:
    samples_with_lane_and_run = f.read().splitlines()

# Python list with sample names only
# This will serve as the list of targets for post-merging steps
with open("processed_data/sample_list.txt") as f:
    samples = f.read().splitlines()

######################
## HELPER FUNCTIONS ##
######################

# A function to create a readgroup for BWA from the sample, lane, and run info
def make_RG(wildcards):
    # Extract sample ID, lane, and run (seq1,seq2, seq3) from input
    sample = wildcards.sample_with_extra.split("_")[0]
    lane = wildcards.sample_with_extra.split("_")[1]
    run = wildcards.sample_with_extra.split("_")[2].replace("rec", "")

    # Assemble the RG header. Fields:
    # ID: ID, the run
    # LB: library, sample + "lib1"
    # PL: platform, ILLUMINA
    # SM: sample, sample
    # PU: platform unit, run + lane + sample

    rg_out = r"@RG\tID:" + run + r"\tLB:" + sample + r"_" + r"lib1" + r"\tPL:ILLUMINA" + r"\tSM:" + sample + r"\tPU:" + run + r"_" + lane + r"_" + sample + r"_lib1" 
  
    return rg_out

# A function to get the input filenames for sample merging
def sample_merge_inputs(sample, directory):
    bases = [x for x in samples_with_lane_and_run if sample in x]
    return [directory + x + "_sorted.bam" for x in bases]


####################
## PIPELINE START ##
####################

rule all:
    input:
        # Multiqc report of trimming data
        "results_prev_seq/multiqc/multiqc_report_trimming.html",
        # MultiQC report of the raw bams
        "results_prev_seq/multiqc/raw_bams/multiqc_report_realigned_sorted_bams.html",
        # Muiltiqc report of the Picard remove dups logs
        "results_prev_seq/multiqc/dedup/multiqc_report_dedup.html",
        # Multiqc report of the final, deduped, realigned, sorted bams
        "results_prev_seq/multiqc/final_bams/multiqc_report_realigned_sorted_bams.html"
        


## trim_raw_reads : remove adaptors and low-quality bases
# Uses fastp. Options:
# -m Merge mode: merge overlapping read pairs which overlap
# -c Correct mismatched bases in the merged region
# Using default parameters for merge and correction
# --detect_adapter_for_pe Auto-detects possible adaptor sequences for removal
# --cut_front Do a sliding window analysis from the 5' end, cut read when quality falls below thresh
# --cut_front_window_size Window size of 5 for sliding window 
# --cut_front_mean_quality Mean quality of 20 for sliding window
# -l 25 Minimum length of 25 BP for the read
# -w Number of cores
# -h, -j Name of report HTML and JSON  output reports
rule trim_and_merge_raw_reads:
    input:
        raw_r1="processed_reads_prev_seq/raw_reads/{sample_with_extra}_R1.fq.gz",
        raw_r2="processed_reads_prev_seq/raw_reads/{sample_with_extra}_R2.fq.gz",
    output:
        trim_merged= "processed_reads_prev_seq/trimmed/{sample_with_extra}.merged.fq.gz",
        trim_r1_pair= "processed_reads_prev_seq/trimmed/{sample_with_extra}.nomerge.pair.R1.fq.gz",
        trim_r2_pair= "processed_reads_prev_seq/trimmed/{sample_with_extra}.nomerge.pair.R2.fq.gz",
        trim_r1_nopair= "processed_reads_prev_seq/trimmed/{sample_with_extra}.nopair.R1.fq.gz",
        trim_r2_nopair= "processed_reads_prev_seq/trimmed/{sample_with_extra}.nopair.R2.fq.gz",
        rep_html= "logs/prev_seq/fastp/{sample_with_extra}_trim_fastp.html",
        rep_json= "logs/prev_seq/fastp/{sample_with_extra}_trim_fastp.json"
    resources:
        cpus = 16
    log:
        "logs/prev_seq/fastp/{sample_with_extra}_trim_log.txt"
    shell:
        """
        fastp -i {input.raw_r1} -I {input.raw_r2} -m --merged_out {output.trim_merged} --out1 {output.trim_r1_pair} --out2 {output.trim_r2_pair} --unpaired1 {output.trim_r1_nopair} --unpaired2 {output.trim_r2_nopair} --detect_adapter_for_pe --cut_front --cut_front_window_size 5 --cut_front_mean_quality 20 -l 25 -j {output.rep_json} -h {output.rep_html} -w $SLURM_CPUS_PER_TASK
        """


## multiqc_trim_reports: collate fastp trimming reports
rule multiqc_trim_reports:
    input:
        expand("logs/prev_seq/fastp/{sample_with_extra}_trim_fastp.json", sample_with_extra = samples_with_lane_and_run)
    output:
        "results_prev_seq/multiqc/multiqc_report_trimming.html"
    params:
        dir_in = "logs/prev_seq/fastp",
        dir_out = "results_prev_seq/multiqc"
    log:
        "logs/prev_seq/multiqc/multiqc_trim_reports.log"
    shell:
        """
        multiqc -f {params.dir_in} -o {params.dir_out} -n multiqc_report_trimming.html > {log} 2>&1
        """

## map_merged_reads: map trimmed, merged reads to reference
#   BWA mem algorithm. Settings:
#   -M Mark shorter split hits as secondary (for Picard compatibility).
#   -t number of threads
#   -R read group, added through lambda function
#   then uses samtools view and samtools sort to convert to bam and sort
#   samtools view options:
#   -b output in bam format
rule map_merged_reads:
    input:
        reads="processed_reads_prev_seq/trimmed/{sample_with_extra}.merged.fq.gz",
        genome=REF,
        genome_index=index_path
    output:
        "processed_reads_prev_seq/mapped/{sample_with_extra}.merged.sorted.bam"
    params:
        basename="processed_reads_prev_seq/mapped/{sample_with_extra}",
        read_group=make_RG
    log:
        "logs/prev_seq/mapping/{sample_with_extra}_merged.log"
    resources:
        cpus=12
    shell:
        """
        # Run bwa mem, pipe to samtools view to convert to bam, pipe to samtools sort
        bwa mem -M -t $SLURM_CPUS_PER_TASK -R '{params.read_group}' {input.genome} {input.reads} 2> {log} | samtools view -b - 2>> {log} | samtools sort - -o {output} 2>> {log}
        """

# # map_unmerged_pairs: map trimmed, not merged, paired reads to reference
#   BWA mem algorithm. Settings:
#   -M Mark shorter split hits as secondary (for Picard compatibility).
#   -t number of threads
#   -R read group, added through lambda function
#   then uses samtools view and samtools sort to convert to bam and sort
#   samtools view options:
#   -b output in bam format
rule map_unmerged_pairs:
    input:
        reads_forward="processed_reads_prev_seq/trimmed/{sample_with_extra}.nomerge.pair.R1.fq.gz",
        reads_reverse="processed_reads_prev_seq/trimmed/{sample_with_extra}.nomerge.pair.R2.fq.gz",
        genome=REF,
        genome_index=index_path
    output:
        "processed_reads_prev_seq/mapped/{sample_with_extra}.nomerge.paired.sorted.bam"
    params:
        basename="processed_reads_prev_seq/mapped/{sample_with_extra}",
        read_group=make_RG
    log:
        "logs/prev_seq/mapping/{sample_with_extra}_nomerge_paired.log"
    resources:
        cpus=12
    shell:
        """
        # Run bwa mem, pipe to samtools view to convert to bam, pipe to samtools sort 
        bwa mem -M -t $SLURM_CPUS_PER_TASK -R '{params.read_group}' {input.genome} {input.reads_forward} {input.reads_reverse} 2> {log} | samtools view -b - 2>> {log} | samtools sort - -o {output} 2>> {log}
        """


## map_unmerged_unpaired: map trimmed, unmerged, unpaired reads to reference
#   BWA mem algorithm. Settings:
#   -M Mark shorter split hits as secondary (for Picard compatibility).
#   -t number of threads
#   -R read group, added through lambda function
#   then uses samtools view and samtools sort to convert to bam and sort
#   samtools view options:
#   -b output in bam format
rule map_unmerged_unpaired:
    input:
        reads_forward="processed_reads_prev_seq/trimmed/{sample_with_extra}.nopair.R1.fq.gz",
        reads_reverse="processed_reads_prev_seq/trimmed/{sample_with_extra}.nopair.R2.fq.gz",
        genome=REF,
        genome_index=index_path
    output:
        mapped_forward = "processed_reads_prev_seq/mapped/{sample_with_extra}.nopair.R1.sorted.bam",
        mapped_reverse = "processed_reads_prev_seq/mapped/{sample_with_extra}.nopair.R2.sorted.bam"
    params:
        basename="processed_reads_prev_seq/mapped/{sample_with_extra}",
        read_group=make_RG
    log:
        forward="logs/prev_seq/mapping/{sample_with_extra}_nopair_R1.log",
        rev="logs/prev_seq/mapping/{sample_with_extra}_nopair_R2.log"
    resources:
        cpus=12
    shell:
        """
        # Run bwa mem, pipe to samtools view to convert to bam, save as a tmp.bam
        # Read 1
        bwa mem -M -t $SLURM_CPUS_PER_TASK -R '{params.read_group}' {input.genome} {input.reads_forward} 2> {log.forward} | samtools view -b - 2>> {log.forward} | samtools sort - -o {output.mapped_forward} 2>> {log.forward}

        # Read 2
        bwa mem -M -t $SLURM_CPUS_PER_TASK -R '{params.read_group}' {input.genome} {input.reads_reverse} 2> {log.rev} | samtools view -b - 2>> {log.rev} | samtools sort - -o {output.mapped_reverse} 2>> {log.rev}
        """

## merge_per_run_bams : merge bam files by sample and run
# merges bams across the 4 types of mapped reads (assembled, paired unassembled, and unpaired SEs)
# for a given sample/lane/sequencing run combination
# use samtools merge, -t is threads, rest is default
rule merge_per_run_bams:
    input: 
        merged="processed_reads_prev_seq/mapped/{sample_with_extra}.merged.sorted.bam",
        unmerge_pair="processed_reads_prev_seq/mapped/{sample_with_extra}.nomerge.paired.sorted.bam",
        unpair_fwd="processed_reads_prev_seq/mapped/{sample_with_extra}.nopair.R1.sorted.bam",
        unpair_rev="processed_reads_prev_seq/mapped/{sample_with_extra}.nopair.R2.sorted.bam"
    log:
        "logs/prev_seq/merge_by_run/{sample_with_extra}_merge.log"
    resources:
        cpus=24
    output:
        temporary("processed_reads_prev_seq/merged_by_run/{sample_with_extra}_sorted.bam")
    shell:
        """
        samtools merge -t {resources.cpus} {output} {input.merged} {input.unmerge_pair} {input.unpair_fwd} {input.unpair_rev} 2> {log}
        """
        
## merge_by_sample: merge all bam files for a given sample
# merges across sequencing runs and lanes.
# file numbers and names vary, are generated by the sample_merge_inputs function
rule merge_by_sample:
    input:
        lambda wildcards: sample_merge_inputs(wildcards.sample, "processed_reads_prev_seq/merged_by_run/")
    output:
        "processed_reads_prev_seq/merged_by_sample/{sample}_sorted.bam"
    log:
        "logs/prev_seq/merge_by_sample/{sample}_merge.log"
    resources:
        cpus=24
    shell:
        """
        samtools merge -t {resources.cpus} {output} {input} 2> {log}
        """

## qualimap_raw_bam: run qualimap on final bam file
# default options, only changed number of threads with -nt
rule qualimap_raw_bam:
    input:
        bam="processed_reads_prev_seq/merged_by_sample/{sample}_sorted.bam"
    output:
        "processed_reads_prev_seq/qc/qualimap_raw/{sample}/qualimapReport.html"
    resources:
        cpus=24
    shell:
        """
        qualimap bamqc -bam {input.bam} -nt {resources.cpus} -outdir processed_reads_prev_seq/qc/qualimap_raw/{wildcards.sample}/ -outformat html --java-mem-size=8G
        """

rule multiqc_raw_bam_reports:
    input:
        expand("processed_reads_prev_seq/qc/qualimap_raw/{sample}/qualimapReport.html",
               sample = samples)
    output:
        "results_prev_seq/multiqc/raw_bams/multiqc_report_realigned_sorted_bams.html"
    shell:
        """
        multiqc processed_reads_prev_seq/qc/qualimap_raw -o results_prev_seq/multiqc/raw_bams -n multiqc_report_realigned_sorted_bams.html
        """

## remove_duplicates: remove duplicates with Picard
# Remove duplicates with Picard. options:
# REMOVE_DUPLICATES=true; remove duplicates, instead of default behavior (which outputs them and flags duplicated reads)
rule remove_duplicates:
    input:
        "processed_reads_prev_seq/merged_by_sample/{sample}_sorted.bam"
    output:
        bam="processed_reads_prev_seq/dedup/{sample}_sorted.bam",
        metrics="logs/prev_seq/remove_duplicates/{sample}_dedup_metrics.txt"
    log:
        log="logs/prev_seq/remove_duplicates/{sample}_dedup_log.log"
    resources:
        cpus=1
    shell:
        """
        picard MarkDuplicates INPUT={input} METRICS_FILE={output.metrics} TMP_DIR=/mnt/beegfs/tt164677e/tempdir OUTPUT={output.bam} REMOVE_DUPLICATES=true 2> {log.log}
        """

## index_deduped_bams: create bam indices for later realignment
rule index_dedup_bams:
    input:
        "processed_reads_prev_seq/dedup/{sample}_sorted.bam"
    output:
        "processed_reads_prev_seq/dedup/{sample}_sorted.bam.bai"
    resources:
        cpus=1
    shell:
        """
        samtools index -b {input}
        """

## multiqc_dedup_reports: collate picard reports of dedup process
rule multiqc_dedup_reports:
    input:
        expand("logs/prev_seq/remove_duplicates/{sample}_dedup_metrics.txt",
        sample = samples)
    output:
        "results_prev_seq/multiqc/dedup/multiqc_report_dedup.html"
    shell:
        """
        multiqc logs/prev_seq/remove_duplicates -o results_prev_seq/multiqc/dedup -n multiqc_report_dedup.html
        """

## realign_target_create: create realignment intervals with GATK3
# default options, only changing number of threads with -nt 
rule realign_target_create:
    input:
        bam="processed_reads_prev_seq/dedup/{sample}_sorted.bam",
        index="processed_reads_prev_seq/dedup/{sample}_sorted.bam.bai",
        genome=REF
    output:
        intervals="processed_reads_prev_seq/realigned/intervals/{sample}_realigner.intervals"
    log:
        create_intervals="logs/prev_seq/realign/{sample}_create_intervals.log"
    resources:
        cpus=24
    shell:
        """
        # Create intervals
        gatk -T RealignerTargetCreator -R {input.genome} -I {input.bam} -o {output.intervals} -nt {resources.cpus} 2> {log.create_intervals}
        """

## realign_indels: realign indels with GATK3
# default options, can't parallelize
rule realign_indels:
    input:
        bam="processed_reads_prev_seq/dedup/{sample}_sorted.bam",
        index="processed_reads_prev_seq/dedup/{sample}_sorted.bam.bai",
        intervals="processed_reads_prev_seq/realigned/intervals/{sample}_realigner.intervals",
        genome=REF
    output:
        "processed_reads_prev_seq/realigned/{sample}_realigned_sorted.bam"
    log:
        "logs/prev_seq/realign/{sample}_realign.log"
    shell:
        """
        gatk -T IndelRealigner -R {input.genome} -I {input.bam} -targetIntervals {input.intervals} -o {output} 2> {log}
        """

## qualimap_final_bam: run qualimap on final bam file
# default options, only changed number of threads with -nt
rule qualimap_final_bam:
    input:
        bam="processed_reads_prev_seq/realigned/{sample}_realigned_sorted.bam"
    output:
        "processed_reads_prev_seq/qc/qualimap_realigned/{sample}/qualimapReport.html"
    resources:
        cpus=24
    shell:
        """
        qualimap bamqc -bam {input.bam} -nt {resources.cpus} -outdir processed_reads_prev_seq/qc/qualimap_realigned/{wildcards.sample}/ -outformat html --java-mem-size=8G
        """



## multiqc_final_bam_reports: collate qualimap reports on final bams
rule multiqc_final_bam_reports:
    input:
        expand("processed_reads_prev_seq/qc/qualimap_realigned/{sample}/qualimapReport.html",
               sample = samples)
    output:
        "results_prev_seq/multiqc/final_bams/multiqc_report_realigned_sorted_bams.html"
    shell:
        """
        multiqc processed_reads_prev_seq/qc -o results_prev_seq/multiqc/final_bams -n multiqc_report_realigned_sorted_bams.html
        """
        
# --- Misc --- #
## dag               : makes a DAG graph for this pipeline
rule dag:
    input: "sf_process_reads_pipeline.smk"
    output: "DAG_process_reads.png"
    shell:
        "snakemake --dag -s sf_process_reads_pipeline.smk --profile sp_genomics | dot -Tpng > DAG_process_reads.png"

## help               : prints help comments for this snakefile
rule help:
    input: "sf_process_reads_pipeline.smk"
    shell:
        "sed -n 's/^##//p' {input}"
