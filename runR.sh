#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 120G
#SBATCH --time=24:00:00
#SBATCH --partition=mrcieu
set -euo pipefail

module load languages/r/3.6.0

# Run R script
Rscript $@