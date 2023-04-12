#!/bin/bash -l
#SBATCH --job-name=single_atac    # Job name
#SBATCH --mem=128G                     # Job memory request
#SBATCH --cpus-per-task=12
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=single_atac.log   # Standard output and error log

module load R
Rscript find_atac_dar.R

