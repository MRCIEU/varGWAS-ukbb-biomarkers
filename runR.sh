#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 24G
#SBATCH --time=24:00:00
#SBATCH --partition=mrcieu,cpu
set -euo pipefail

module load languages/r/3.6.0
module load apps/plink/1.90
PATH=$PATH:/mnt/storage/home/ml18692/scratch/db/locuszoom/locuszoom/bin

# Run R script
Rscript $@