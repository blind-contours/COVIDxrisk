#!/bin/bash
# Job name:
#SBATCH --job-name=day100_cases_COVID
#
# Partition:
#SBATCH --partition=savio3_bigmem
#
#SBATCH --qos=savio_lowprio
#SBATCH --account=co_biostat
#
# Number of nodes for use case:
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --exclusive
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=0
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=david_mccoy@berkeley.edu
## Command(s) to run (example):
module load r/4.0.3
cd /global/scratch/users/david_mccoy/COVIDxrisk/COVIDxrisk
R CMD BATCH --no-save \
  R/03_fit_day100_cases_model.R \
  logs/cases_day100_modeling.Rout
