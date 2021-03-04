#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 48G
#SBATCH --time=120:00:00
#SBATCH --partition=mrcieu
set -euo pipefail

module load build/gcc-5.5.0

# set logging level
export SPDLOG_LEVEL=debug

# run vGWAS using B-P
/mnt/storage/home/ml18692/projects/jlst_cpp/build/bin/jlst_cpp \
-v data/ukb_bmi.txt \
-s , \
-o data/ukb_bmi.vgwas.t1.txt \
-b /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen \
-p body_mass_index.21001.0.0 \
-c sex.31.0.0,age_at_recruitment.21022.0.0,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
-i appieu \
-t 4
