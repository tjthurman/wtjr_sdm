#!/bin/bash

#SBATCH --nodes=1 # For programs that aren't running with MPI parallelization
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # the number of cores to use. If a serial job, do 1.
#SBATCH --output="/home/tt164677e/tim_beegfs/wtjr_genomics/logs/%x_%j.out" #where to send the errors and output from the console.
#SBATCH --mail-user=timothy.j.thurman@gmail.com #where to send email to user
#SBATCH --mail-type=ALL #when to email the user
#SBATCH --partition=good_lab_cpu
#SBATCH --job-name=realSFS_fst_print

# Run with Fst conda environment

cd /home/tt164677e/tim_beegfs/wtjr_genomics/results/angsd_Fst/colorado_samples_bamlist_fixed_NDK_samples_bamlist/wtjr_COL-NDK_20221113/by_scaff_idx


for file in *.fst.idx
do
realSFS fst print $file >> ../genome_fst_print.txt
done 