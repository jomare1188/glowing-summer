#!/bin/sh
#SBATCH --job-name=qc_illumina
#SBATCH --partition=long
#SBATCH --ntasks-per-node=30
#SBATCH --mem=50gb
#SBATCH --error=%j_qc.err
#SBATCH --output=%j_qc.out


source ~/.bashrc
conda activate qc_illumina

bash qc_pipeline.sh /home/dmpachon/jorge/TATIANA/fastq "" 30



