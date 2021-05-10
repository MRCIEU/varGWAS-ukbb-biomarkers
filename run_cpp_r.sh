#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem 100G
#SBATCH --time=100:00:00
#SBATCH --partition=mrcieu
set -euo pipefail

module load build/gcc-5.5.0

# set trait from args
trait="$1"

# set chr from args
chr="$2"

# run vGWAS
/mnt/storage/home/ml18692/projects/jlst_cpp/build/bin/jlst_cpp \
-v data/"$trait".txt \
-s , \
-o data/"$trait".vgwas.chr"$chr".txt \
-b /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr"$chr".bgen \
-p "$trait" \
-c sex.31.0.0,age_at_recruitment.21022.0.0,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
-i appieu \
-r \
-t 27