#!/bin/bash
#SBATCH --job-name=run_12m
#SBATCH --output=run12m_%A_%a.out
#SBATCH --error=run12m_%A_%a.err
#SBATCH --time=168:00:00
#SBATCH --mem=400G
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=12


# Load modules
module load R

#/usr/bin/time -v Rscript tmp_old_p01m.R
#/usr/bin/time -v Rscript tmp_new_p01m.R
#/usr/bin/time -v Rscript tmp_old_p1m.R
#/usr/bin/time -v Rscript tmp_old_p1m.R
/usr/bin/time -v Rscript tmp_new_m.R

