#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 120G
#SBATCH --time=24:00:00
#SBATCH --partition=mrcieu
set -euo pipefail

module load languages/anaconda3/2020.02-tflow-2.2.0
source activate ldsc
PATH=$PATH:/mnt/storage/home/ml18692/scratch/apps/ldsc

# munge GWAS
munge_sumstats.py \
--sumstats "alanine_aminotransferase.30620.0.0.vgwas.chr01.txt" \
--out "alanine_aminotransferase.30620.0.0" \



ldsc.py -h