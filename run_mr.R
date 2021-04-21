#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 8G
#SBATCH --time=20:00:00
#SBATCH --partition=mrcieu
set -euo pipefail

module load languages/r/3.6.0

# find disease outcomes associated with vQTL
Rscript mr.R \
-t "$1" \
-i "$2" \
-o data/"$1".mr.txt