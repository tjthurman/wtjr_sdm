# A pipeline to process raw read data
# of WTJR across North America into mapped, filtered .bam files

# This pipeline processes the files from my (TJT) sequencing
# There's a separate file for processing Mafalda's, using similar settings but requiring
# some mergine steps, as some individuals were sequenced across multiple lanes


# Python env setup
import os
import re
# import csv
# import sys
# from datetime import datetime

###############
##  GLOBALS  ##
###############

# Reference genome
REF="/mnt/beegfs/tt164677e/genomes/lepus_townsendii/DMNS18807_06042020_pseudohap2.1_10k.fasta"


# Folder with raw read data for my sequencing
RAW_DATA_FOLDER = "/mnt/beegfs/tt164677e/wtjr_genomics/raw_read_data/usftp21.novogene.com/raw_data"

# Path to genome index file for BWA
index_path= REF + ".amb"

###############
##   SETUP   ##
############### 

# Pull all sample files from the raw data folder
# get filenames of the individual fastas
# removed the undetermined reads for now, could come back to try to recover those

fasta_fullpath = []
samples = []
for root, dirs, files in os.walk(RAW_DATA_FOLDER):
    for name in files:
        if re.search("_CKDL\d+-1a-\S+_HH5YYDSX2_L1_[12].fq.gz$", name):
            samples.append(re.sub("_CKDL\d+-1a-\S+_HH5YYDSX2_L1_[12].fq.gz$", "", name))
            fasta_fullpath.append(os.path.join(root, name))
            
# Get unique sample IDs
samples = list(set(samples))
samples.sort()

# Subset down to just a few for now, for testing
# samples = samples[1:5]



######################
## HELPER FUNCTIONS ##
######################

# At the very start, need to match each
# sample ID back to its fastq files.
# This assumes that each sample only has one set of sequence files in a given folder. 
def get_R1_for_sample(wildcards):
    outfile = "file_not_found.txt"
    for filename in fasta_fullpath:
        if re.search(wildcards.sample, filename):
            if re.search("L1_1", filename):
                outfile = filename
    return outfile

def get_R2_for_sample(wildcards):
    outfile = "file_not_found.txt"
    for filename in fasta_fullpath:
        if re.search(wildcards.sample, filename):
            if re.search("L1_2", filename):
                outfile = filename
    return outfile

def make_RG(wildcards):
    for filename in fasta_fullpath:
        if re.search(wildcards.sample, filename):
            if re.search("L1_1", filename):
                basename = os.path.basename(filename)
    # Extract sample ID, lane, and run (seq1,seq2, seq3) from input
    sample_ID = re.sub("_CKDL\d+", "", basename.split("-")[0])
    library = re.sub("_HH5", "", re.findall("CKDL\S+_HH5", basename)[0])
    flowcell = basename.split("-")[3].split("_")[1]
    lane = re.findall("L\d", basename.split("-")[3])[0]
    # Assemble the RG header. Fields:
    # ID: Individual sample ID plus sample number
    # LB: library, sample + "lib1"
    # PL: platform, ILLUMINA
    # SM: sample, sample
    # PU: platform unit, run + lane + sample
    rg_out = "@RG\\tID:" + sample_ID  + "\\tLB:" + library + "\\tPL:ILLUMINA" + "\\tSM:" + sample_ID + "\\tPU:" + "X202SC21062583" + "-" + lane + "-" + sample_ID
    return rg_out

####################
## PIPELINE START ##
####################
localrules: all

# Onstart and onsuccess, make a log files and copy the snakefile that was used
# into the results directory for posterity
# pipeline_log_file = bd("logs/pipeline_log.tsv")
# onstart:
#     os.makedirs(os.path.dirname(pipeline_log_file), exist_ok=True)
#     with open(pipeline_log_file, 'w') as tsvfile:
#         writer = csv.writer(tsvfile, delimiter='\t')
#         writer.writerow(["Raw data folder used:", RAW_DATA_FOLDER])
#         writer.writerow(["Reference genome used:", REF])
#         writer.writerow(["GFF file used:", GFF])
#         writer.writerow(["Start time:", start_time.strftime("%B %d, %Y: %H:%M:%S")])
        

# onsuccess:
#     smk_copy_command = 'cp pipe_reads_to_lineages.smk ' + str(bd("snakefile_used_copy.smk"))
#     end_time = datetime.now()
#     elapsed_time = end_time - start_time
#     os.popen(smk_copy_command) 
#     with open(pipeline_log_file, 'a') as tsvfile:
#         writer = csv.writer(tsvfile, delimiter='\t')
#         writer.writerow(["End time:", end_time.strftime("%B %d, %Y: %H:%M:%S")])
#         writer.writerow(["Elapsed time:", str(elapsed_time)])
#         writer.writerow(["Copy of snakefile used stored at:", str(bd("snakefile_used_copy.smk"))])

# all: The rule that looks for the final desired output files to initiate running all rules to generate those files.
rule all:
    input:
        "results_new_seq/multiqc/multiqc_report_trimming.html",
        "results_new_seq/multiqc/multiqc_report_raw_bams.html",
        "results_new_seq/multiqc/multiqc_report_dedup.html",
        "results_new_seq/multiqc/multiqc_report_dedup_realigned_sorted_bams.html"



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
        raw_r1=get_R1_for_sample,
        raw_r2=get_R2_for_sample
    output:
        trim_merged= "processed_reads_new_seq/trimmed/{sample}.merged.fq.gz",
        trim_r1_pair= "processed_reads_new_seq/trimmed/{sample}.nomerge.pair.R1.fq.gz",
        trim_r2_pair= "processed_reads_new_seq/trimmed/{sample}.nomerge.pair.R2.fq.gz",
        trim_r1_nopair= "processed_reads_new_seq/trimmed/{sample}.nopair.R1.fq.gz",
        trim_r2_nopair= "processed_reads_new_seq/trimmed/{sample}.nopair.R2.fq.gz",
        rep_html= "logs/new_seq/fastp/{sample}_trim_fastp.html",
        rep_json= "logs/new_seq/fastp/{sample}_trim_fastp.json"
    resources:
        cpus = 16
    log:
        "logs/new_seq/fastp/{sample}_trim_log.txt"
    shell:
        """
        fastp -i {input.raw_r1} -I {input.raw_r2} -m --merged_out {output.trim_merged} --out1 {output.trim_r1_pair} --out2 {output.trim_r2_pair} --unpaired1 {output.trim_r1_nopair} --unpaired2 {output.trim_r2_nopair} --detect_adapter_for_pe --cut_front --cut_front_window_size 5 --cut_front_mean_quality 20 -l 25 -j {output.rep_json} -h {output.rep_html} -w $SLURM_CPUS_PER_TASK 2> {log}
        """


## multiqc_trim_reports: collate fastp trimming reports
rule multiqc_trim_reports:
    input:
        expand("logs/new_seq/fastp/{sample}_trim_fastp.json", sample = samples)
    output:
        "results_new_seq/multiqc/multiqc_report_trimming.html"
    params:
        dir_in = "logs/new_seq/fastp",
        dir_out = "results_new_seq/multiqc"
    log:
        "logs/new_seq/multiqc/multiqc_trim_reports.log"
    shell:
        """
        multiqc -f {params.dir_in} -o {params.dir_out} -n multiqc_report_trimming.html > {log} 2>&1
        """

## index_ref: index genome for BWA
# Index the reference genome, if it isn't already
rule index_ref:
    input:
        REF
    output:
        multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/new_seq/index_ref.log"
    shell:
        """
        bwa index {input} > {log} 2>&1
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
        reads="processed_reads_new_seq/trimmed/{sample}.merged.fq.gz",
        genome=REF,
        genome_index=index_path
    output:
        "processed_reads_new_seq/mapped/{sample}.merged.sorted.bam"
    params:
        basename="processed_reads_new_seq/mapped/{sample}",
        read_group=make_RG
    log:
        "logs/new_seq/mapping/{sample}_merged.log"
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
        reads_forward="processed_reads_new_seq/trimmed/{sample}.nomerge.pair.R1.fq.gz",
        reads_reverse="processed_reads_new_seq/trimmed/{sample}.nomerge.pair.R2.fq.gz",
        genome=REF,
        genome_index=index_path
    output:
        "processed_reads_new_seq/mapped/{sample}.nomerge.paired.sorted.bam"
    params:
        basename="processed_reads_new_seq/mapped/{sample}",
        read_group=make_RG
    log:
        "logs/new_seq/mapping/{sample}_nomerge_paired.log"
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
        reads_forward="processed_reads_new_seq/trimmed/{sample}.nopair.R1.fq.gz",
        reads_reverse="processed_reads_new_seq/trimmed/{sample}.nopair.R2.fq.gz",
        genome=REF,
        genome_index=index_path
    output:
        mapped_forward = "processed_reads_new_seq/mapped/{sample}.nopair.R1.sorted.bam",
        mapped_reverse = "processed_reads_new_seq/mapped/{sample}.nopair.R2.sorted.bam"
    params:
        basename="processed_reads_new_seq/mapped/{sample}",
        read_group=make_RG
    log:
        forward="logs/new_seq/mapping/{sample}_nopair_R1.log",
        rev="logs/new_seq/mapping/{sample}_nopair_R2.log"
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

## merge_bams_by_sample : merge bam files by sample and run
# merges bams across the 4 types of mapped reads (assembled, paired unassembled, and unpaired SEs)
# for a given sample/lane/sequencing run combination
# use samtools merge, -t is threads
# -c merges identical readgroup headers, which our files from the same individual should have. 
rule merge_sample_bams:
    input: 
        merged="processed_reads_new_seq/mapped/{sample}.merged.sorted.bam",
        unmerged_pair="processed_reads_new_seq/mapped/{sample}.nomerge.paired.sorted.bam",
        nopair_fwd="processed_reads_new_seq/mapped/{sample}.nopair.R1.sorted.bam",
        nopair_rev="processed_reads_new_seq/mapped/{sample}.nopair.R2.sorted.bam"
    log:
        "logs/new_seq/merge_bams/{sample}_merge.log"
    resources:
        cpus=4
    output:
        "processed_reads_new_seq/per_sample_bams/{sample}.sorted.bam"
    shell:
        """
        samtools merge -c -t {resources.cpus} {output} {input.merged} {input.unmerged_pair} {input.nopair_fwd} {input.nopair_rev} 2> {log}
        """

## index_raw_bams: index bams
rule index_raw_bams:
    input:
        "processed_reads_new_seq/per_sample_bams/{sample}.sorted.bam"
    output:
        "processed_reads_new_seq/per_sample_bams/{sample}.sorted.bam.bai"
    log:
        "logs/new_seq/index_bams/{sample}.log"
    shell:
        """
        samtools index -b {input} 2> {log}
        """

## qualimap_raw_bam: run qualimap on raw bam file
# default options, only changed number of threads with -nt
rule qualimap_raw_bam:
    input:
        bam="processed_reads_new_seq/per_sample_bams/{sample}.sorted.bam",
        bai="processed_reads_new_seq/per_sample_bams/{sample}.sorted.bam.bai"
    output:
        "processed_reads_new_seq/QC/qualimap/{sample}/qualimapReport.html"
    params:
        out_dir="processed_reads_new_seq/QC/qualimap/{sample}"
    resources:
        cpus=8
    log:
        "logs/new_seq/qualimap/{sample}.log"
    shell:
        """
        qualimap bamqc -bam {input.bam} -nt $SLURM_CPUS_PER_TASK -outdir {params.out_dir} -outformat html --java-mem-size=4G > {log} 2>&1
        """

## multiqc_raw_bam_report: collate qualimap reports on raw bams
rule multiqc_raw_bam_report:
    input:
        expand("processed_reads_new_seq/QC/qualimap/{sample}/qualimapReport.html", sample = samples)
    output:
        "results_new_seq/multiqc/multiqc_report_raw_bams.html",
        "results_new_seq/multiqc/multiqc_report_raw_bams_data/multiqc_general_stats.txt"
    params:
        dir_in = "processed_reads_new_seq/QC/qualimap",
        dir_out = "results_new_seq/multiqc"
    log:
        "logs/new_seq/multiqc/multiqc_raw_bam_reports.log"
    shell:
        """
        multiqc -f {params.dir_in} -o {params.dir_out} -n multiqc_report_raw_bams.html > {log} 2>&1
        """


## remove_duplicates: remove duplicates with Picard
# Remove duplicates with Picard. options:
# REMOVE_DUPLICATES=true; remove duplicates, instead of default behavior (which outputs them and flags duplicated reads)
rule remove_duplicates:
    input:
        bam="processed_reads_new_seq/per_sample_bams/{sample}.sorted.bam",
        bai="processed_reads_new_seq/per_sample_bams/{sample}.sorted.bam.bai"
    output:
        bam="processed_reads_new_seq/dedup/{sample}.sorted.bam",
        metrics="logs/new_seq/remove_duplicates/{sample}_dedup_metrics.txt"
    log:
        log="logs/new_seq/remove_duplicates/{sample}_dedup_log.log"
    resources:
        cpus=1
    shell:
        """
        picard MarkDuplicates INPUT={input.bam} METRICS_FILE={output.metrics} OUTPUT={output.bam} REMOVE_DUPLICATES=true 2> {log.log}
        """

## index_deduped_bams: create bam indices for later realignment
rule index_dedup_bams:
    input:
        "processed_reads_new_seq/dedup/{sample}.sorted.bam"
    output:
        "processed_reads_new_seq/dedup/{sample}.sorted.bam.bai"
    resources:
        cpus=1
    shell:
        """
        samtools index -b {input}
        """

## multiqc_dedup_reports: collate picard reports of dedup process
rule multiqc_dedup_reports:
    input:
        expand("logs/new_seq/remove_duplicates/{sample}_dedup_metrics.txt",
        sample = samples)
    output:
        "results_new_seq/multiqc/multiqc_report_dedup.html"
    params:
        dir_in = "logs/new_seq/remove_duplicates",
        dir_out = "results_new_seq/multiqc"
    log:
        log="logs/new_seq/multiqc/multiqc_dedup_reports.log"
    shell:
        """
        multiqc -f {params.dir_in} -o {params.dir_out} -n multiqc_report_dedup.html > {log} 2>&1
        """

## Install GATK 3
rule install_gatk_3:
    output:
        "bin/GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2"
    shell:
        """
        cd bin
        wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2
        gatk-register GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2 
        """

## realign_target_create: create realignment intervals with GATK3
# default options, only changing number of threads with -nt 
rule realign_target_create:
    input:
        bam="processed_reads_new_seq/dedup/{sample}.sorted.bam",
        index="processed_reads_new_seq/dedup/{sample}.sorted.bam.bai",
        genome=REF
    output:
        intervals="processed_reads_new_seq/realign_targets/{sample}_realigner.intervals"
    log:
        create_intervals="logs/new_seq/realign/{sample}_create_intervals.log"
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
        bam="processed_reads_new_seq/dedup/{sample}.sorted.bam",
        index="processed_reads_new_seq/dedup/{sample}.sorted.bam.bai",
        intervals="processed_reads_new_seq/realign_targets/{sample}_realigner.intervals",
        genome=REF
    output:
        "processed_reads_new_seq/realigned/{sample}.realigned.sorted.bam"
    log:
        "logs/new_seq/realign/{sample}_realign.log"
    shell:
        """
        gatk -T IndelRealigner -R {input.genome} -I {input.bam} -targetIntervals {input.intervals} -o {output} 2> {log}
        """

## qualimap_dedup_realign_bam: run qualimap on deduped, realigned bam file
# default options, only changed number of threads with -nt
rule qualimap_dedup_realign_bam:
    input:
        bam="processed_reads_new_seq/realigned/{sample}.realigned.sorted.bam"
    output:
        "processed_reads_new_seq/QC/qualimap_realigned/{sample}/qualimapReport.html"
    params:
        out_dir="processed_reads_new_seq/QC/qualimap_realigned/{sample}"
    resources:
        cpus=8
    log:
        "logs/new_seq/qualimap_realigned/{sample}.log"
    shell:
        """
        qualimap bamqc -bam {input.bam} -nt $SLURM_CPUS_PER_TASK -outdir {params.out_dir} -outformat html --java-mem-size=4G > {log} 2>&1
        """

## multiqc_final_bam_reports: collate qualimap reports on final bams
rule multiqc_dedup_realigned_bam_report:
    input:
        expand("processed_reads_new_seq/QC/qualimap_realigned/{sample}/qualimapReport.html",
               sample = samples)
    output:
        "results_new_seq/multiqc/multiqc_report_dedup_realigned_sorted_bams.html"
    params:
        dir_in = "processed_reads_new_seq/QC/qualimap_realigned",
        dir_out = "results_new_seq/multiqc"
    log:
        "logs/new_seq/multiqc/multiqc_dedup_realigned_bam_report.log"
    shell:
        """
        multiqc -f {params.dir_in} -o {params.dir_out} -n multiqc_report_dedup_realigned_sorted_bams.html > {log} 2>&1
        """
