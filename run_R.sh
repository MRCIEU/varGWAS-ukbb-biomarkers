#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 32G
#SBATCH --time=200:00:00
#SBATCH --partition=mrcieu
set -euo pipefail

module load languages/r/3.6.0

# run sim
Rscript sim.R --trait alanine_aminotransferase.30620