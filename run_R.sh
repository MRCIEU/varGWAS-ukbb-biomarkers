#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 32G
#SBATCH --time=200:00:00
#SBATCH --partition=mrcieu
set -euo pipefail

module load languages/r/3.6.0

# run vGWAS using B-P
Rscript bp.R \
-p data/ukb_bmi.txt \
-t body_mass_index.21001 \
-s data/snps.txt \
-o data/ukb_bmi.vgwas.r_subsample.txt