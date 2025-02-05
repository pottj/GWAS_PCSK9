#!/bin/bash

## Example SLURM script for BSU icelake jobs

## Section 1: SLURM Commands

## All SLURM commands must be placed at the start of the file
## Full documentation can be found here: https://slurm.schedmd.com/sbatch.html

## Enter a short name for the job, to be shown in SLURM output
#SBATCH -J regenieStep2

## Enter the wall-clock time limit for your jobs.
## If jobs reach this limit they are automatically killed.
## Maximum value 36:00:00.
#SBATCH --time=06:00:00

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
#SBATCH --array=1-22

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

## Step 5b: run regenie step 2 per chromosome
#
# For quantitative traits, we use a linear regression model for association testing.
#
# - Covariates are regressed out of the phenotypes and genetic markers.
# - The LOCO predictions from Step 1 are removed from the phenotypes.
# - Linear regression is then used to test association of the residualized phenotype and the genetic marker.
# - Parallel linear algebra operations in the Eigen library are used where possible.
#

regenie \
  --step 2 \
  --bgen /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c${SLURM_ARRAY_TASK_ID}_b0_v3.bgen \
  --ref-first \
  --sample /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c${SLURM_ARRAY_TASK_ID}_b0_v3_s487160.sample \
  --phenoFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_phenotypes_EUR_step2_PCSK9.txt \
  --covarFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_covariates_EUR.txt \
  --firth --approx --pThresh 0.01 \
  --pred /rds/user/jp2047/hpc-work/GWAS_PCSK9/regenie/ukb_step1_PCSK9_pred.list \
  --bsize 400 \
  --threads 20 \
  --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/regenie/ukb_step2_PCSK9_chr${SLURM_ARRAY_TASK_ID}
 
