#!/bin/bash
# Job name:
#SBATCH --job-name=covid_cases_1year_modeling
#
# Partition:
#SBATCH --partition=savio2
#

#SBATCH --qos=biostat_savio2_normal
#SBATCH --account=co_biostat

# Number of nodes for use case:
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time 4:00:00
#
## Command(s) to run (example):
module load r/4.0.3
OMP_NUM_THREADS=1
### for foreach+doSNOW ###
cd /global/scratch/users/david_mccoy/COVIDxrisk/COVIDxrisk
 
R CMD BATCH --no-save R/03_fit_year1_cases_model.R logs/cases_1year_modeling.Rout