#!/bin/sh
#SBATCH --job-name=snp_calling
#SBATCH --partition=short
#SBATCH --ntasks-per-node=5
#SBATCH --mem=100gb
#SBATCH --error=%j_snp_claude.err
#SBATCH --output=%j_snp_claude.out


source ~/.bashrc
conda activate SNP_call

#bash snp_calling.sh "/home/dmpachon/jorge/TATIANA/qc/qc_results/trimmed_reads/" "/home/dmpachon/jorge/TATIANA/CallicarpaGenome/car_asm.fa" "" 60

#bash snp_call_claude.sh
#bash snp_paralell.sh
bash snp_optimized.sh
