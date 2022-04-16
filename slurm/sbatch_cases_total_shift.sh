#!/bin/bash
# Job name:
#SBATCH --job-name=covid_cases_total_shift
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
#SBATCH --time 6:00:00
#
## Command(s) to run (example):
module load r/4.0.3
OMP_NUM_THREADS=1
### for foreach+doSNOW ###
cd /global/scratch/users/david_mccoy/COVIDxrisk/COVIDxrisk
 
R CMD BATCH --no-save R/04_joint_shift_total_cases.R logs/cases_total_forward_shift.Rout