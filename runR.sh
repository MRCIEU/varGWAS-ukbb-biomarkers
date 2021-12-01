#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 8G
#SBATCH --time=24:00:00
#SBATCH --partition=short
set -euo pipefail

module load lang/r/4.0.3-bioconductor-gcc
module load apps/plink/1.90

# Run R script
Rscript $@
