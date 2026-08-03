#!/bin/sh
#SBATCH --job-name=snp_calling
#SBATCH --partition=long
#SBATCH --ntasks-per-node=80
#SBATCH --mem=500gb
#SBATCH --error=%j_snp_invariants.err
#SBATCH --output=%j_snp_invariants.out


source ~/.bashrc
conda activate SNP_call

#bash snp_calling.sh "/home/dmpachon/jorge/TATIANA/qc/qc_results/trimmed_reads/" "/home/dmpachon/jorge/TATIANA/CallicarpaGenome/car_asm.fa" "" 60

#bash snp_call_claude.sh
#bash snp_paralell.sh
#bash snp_optimized.sh
bash snp_invariants.sh
