#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 120G
#SBATCH --time=24:00:00
#SBATCH --partition=mrcieu
set -euo pipefail

module load languages/anaconda2/5.3.1.tensorflow-1.12

/mnt/storage/scratch/ml18692/.conda/ldsc/bin/python \
/mnt/storage/home/ml18692/apps/ldsc/ldsc.py \
--h2 "$1" \
--ref-ld-chr /mnt/storage/home/ml18692/db/ldsc/eur_w_ld_chr/ \
--w-ld-chr /mnt/storage/home/ml18692/db/ldsc/eur_w_ld_chr/ \
--out "$1"

rm "$1"