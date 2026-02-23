#!/bin/bash
#SBATCH --job-name=run_2
#SBATCH --output=run2_%A_%a.out
#SBATCH --error=run2_%A_%a.err
#SBATCH --time=168:00:00
#SBATCH --mem=400G
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=2

# Load modules
module load R
Rscript run.R

