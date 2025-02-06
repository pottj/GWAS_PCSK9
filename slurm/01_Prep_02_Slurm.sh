#!/bin/bash

## Example SLURM script for BSU icelake jobs

## Section 1: SLURM Commands

## All SLURM commands must be placed at the start of the file
## Full documentation can be found here: https://slurm.schedmd.com/sbatch.html

## Enter a short name for the job, to be shown in SLURM output
#SBATCH -J 2_regenie1

## Enter the wall-clock time limit for your jobs.
## If jobs reach this limit they are automatically killed.
## Maximum value 36:00:00.
#SBATCH --time=00:30:00

## For single-core jobs, this number should be '1'. 
## If your job has built-in parallelism, eg using OpenMP or 
## R's foreach() and doParallel(), increase this number as desired.
## The maximum value is 76 on icelake; 112 on sapphire
#SBATCH --cpus-per-task=20

## Each task is allocated 3.3G (icelake) or 6.7G (icelake-himem) or 4.6G (sapphire)
## If this is insufficient, uncomment and edit this line.
## Maximum value 256G (icelake/sapphire) or 512G (icelake-himem)
#SBATCH --mem=128G

## The system can send emails when your job starts and stops.
## Values include BEGIN, END, ALL, and TIME_LIMIT_80 and TIME_LIMIT_90 
## (reaching 80% or 90% of time limit.) Specify ARRAY_TASKS to receive
## a separate mail for each task. Multiple values can be given, separated by a comma.
#SBATCH --mail-type=FAIL

## The project account name.
## Use mrc-bsu-sl2-cpu for icelake and mrc-bsu-sl2-gpu for ampere
#SBATCH -A mrc-bsu-sl2-cpu

## The partition. Use icelake for normal jobs, or icelake-himem if needed.
#SBATCH -p icelake

## GPU jobs only:
## Uncomment and specify the number of GPUs required per node, maximum 4.
## Note that there is a maximum of 3 cores per GPU.
## The gpu partition is ampere.
## #SBATCH --gres=gpu:1

## Array jobs:
## Start multiple jobs at once.
## Note that resources (cores, memory, time) requested above are for each
## individual array task, NOT the total array.
#SBATCH --array=1-5

##  - - - - - - - - - - - - - -

## Section 2: Modules

# All scripts should include the first three lines.

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment

# Load the latest R version.
# Before running your code, you should run R and install any required packages.
module load ceuadmin/regenie/3.2.9

# If using the GPU cluster, replace the third line with the uncommented line:
# module load rhel8/default-amp

#! Insert additional module load commands after this line if needed:

## - - - - - - - - - - -

## Section 3: Run your application

## Step 4: run regenie step 1 per phenotype
#
# Whole genome regression model is fit at a subset of the total set of available genetic markers and get Leave One Chromosome Out (LOCO) predictions

regenie \
  --step 1 \
  --bed /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/ukb_cal_allChrs \
  --extract /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/qc_pass.snplist \
  --keep /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/qc_pass.id \
  --phenoFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_phenotypes_EUR_step1_PCSK9_${SLURM_ARRAY_TASK_ID}.txt \
  --covarFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_covariates_EUR.txt \
  --threads 20 \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix /rds/user/jp2047/hpc-work/GWAS_PCSK9/regenie/tmpdir/regenie_tmp_preds \
  --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/regenie/ukb_step1_PCSK9_${SLURM_ARRAY_TASK_ID}
  
