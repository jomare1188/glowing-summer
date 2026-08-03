#!/bin/sh
#SBATCH --job-name=pixy_pi
#SBATCH --partition=long
#SBATCH --ntasks-per-node=50
#SBATCH --mem=100gb
#SBATCH --error=%j_pixy_pi.err
#SBATCH --output=%j_pixy_pi.out


source ~/.bashrc
conda activate SNP_call

#bash snp_calling.sh "/home/dmpachon/jorge/TATIANA/qc/qc_results/trimmed_reads/" "/home/dmpachon/jorge/TATIANA/CallicarpaGenome/car_asm.fa" "" 60

#bash snp_call_claude.sh
#bash snp_paralell.sh
#bash snp_optimized.sh
pixy   --vcf cohort.allsites.vcf.gz   --populations populations.txt   --window_size 10000   --n_cores 50  --output_prefix pixy_out --stats pi --output_folder pixy_out
