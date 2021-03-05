#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 32G
#SBATCH --time=500:00:00
#SBATCH --partition=mrcieu
set -euo pipefail

module load languages/r/3.6.0

# set chr from cmd arg
chr="$1"

# run vGWAS using B-P
Rscript bp.R "$chr"