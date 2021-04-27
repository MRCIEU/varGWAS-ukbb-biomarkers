#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem 8G
#SBATCH --partition=mrcieu
set -euo pipefail

module load languages/r/3.6.0

# set trait from args
trait="$1"

# QC script
Rscript qc.R -t "$trait"