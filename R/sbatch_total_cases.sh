#!/bin/bash
# Job name:
#SBATCH --job-name=total_covid_cases_modeling
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
R CMD BATCH --no-save 03_fit_total_cases_model.R logs/total_cases_modeling.Rout
