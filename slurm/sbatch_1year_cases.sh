#!/bin/bash
# Job name:
#SBATCH --job-name=1year_cases_COVID
#
# Partition:
#SBATCH --partition=savio2
#
#SBATCH --qos=biostat_savio2_normal
#SBATCH --account=co_biostat
#
# Number of nodes for use case:
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --exclusive
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=72:00:00
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
  R/03_fit_year1_cases_model.R \
  logs/cases_1year_modeling.Rout
