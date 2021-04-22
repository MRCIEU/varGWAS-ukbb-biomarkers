#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --partition=veryshort
set -euo pipefail

module load languages/r/3.6.0

# set trait from args
trait="$1"

# set chr from args
chr="$2"

# validate tophits using B-F
Rscript validate.R \
-p data/"$trait".txt \
-t "$trait" \
-g data/"$trait".vgwas.chr"$chr".txt \
-o data/"$trait".validate.chr"$chr".txt