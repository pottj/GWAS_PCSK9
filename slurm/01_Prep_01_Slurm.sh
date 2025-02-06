#!/bin/bash

## Example SLURM script for BSU icelake jobs

## Section 1: SLURM Commands

## All SLURM commands must be placed at the start of the file
## Full documentation can be found here: https://slurm.schedmd.com/sbatch.html

## Enter a short name for the job, to be shown in SLURM output
#SBATCH -J 1_PrepData

## Enter the wall-clock time limit for your jobs.
## If jobs reach this limit they are automatically killed.
## Maximum value 36:00:00.
#SBATCH --time=01:00:00

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
## #SBATCH --array=3-10

##  - - - - - - - - - - - - - -

## Section 2: Modules

# All scripts should include the first three lines.

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment

# Load the latest R version.
# Before running your code, you should run R and install any required packages.
module load R/4.3.1-icelake
module load plink/1.9

# If using the GPU cluster, replace the third line with the uncommented line:
# module load rhel8/default-amp

#! Insert additional module load commands after this line if needed:

## - - - - - - - - - - -

## Section 3: Run your application

## Step 1: create input 
#
# Step 1 and 2 of regenie need the necessary sample files for the --keep, --phenoFile and --covarFile commands

R CMD BATCH --vanilla ../scripts/01_Prep_01_checkUKB_samples.R ../scripts/01_Prep_01_checkUKB_samples.R.out
cp Rplots.pdf 01_Prep_01_Rplots.pdf
rm Rplots.pdf

## Step 2: preparing the genotype file
#
# Step 1 of regenie needs 1 genotype file (not just imputed genedosages)

rm -f /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/list_beds.txt
for chr in {2..22}; do echo "/rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c${chr}_b0_v2.bed /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb_snp_chr${chr}_v2.bim /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2_s488131.fam" >> /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/list_beds.txt; done

plink \
  --bed /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2.bed \
  --bim /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb_snp_chr1_v2.bim \
  --fam /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2_s488131.fam \
  --threads 20 \
  --keep /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_SampleList_EUR.txt \
  --merge-list /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/list_beds.txt \
  --make-bed --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/ukb_cal_allChrs

# Check if exit code != 0 for last command (plink)
retVal =$?
if [ $retVal -ne 0 ]; then
    echo "Error"
    exit $retVal 
fi 

## Step 3: Exclusion files
#
# QC applied to the generated genotype file per SNP (MAF, MAC, missing genotype rates per SNP, and HWE) and samples (missing genotype rates per sample)

module unload plink
module load plink/2.00-alpha

plink2 \
  --bfile /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/ukb_cal_allChrs \
  --keep /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_SampleList_EUR.txt \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --mind 0.1 \
  --threads 20 \
  --write-snplist --write-samples --no-id-header \
  --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/qc_pass

plink2 \
--bfile /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/ukb_cal_allChrs \
--keep /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_SampleList_EUR_extra.txt \
--maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
--mind 0.1 \
--threads 20 \
--write-snplist --write-samples --no-id-header \
--out /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/qc_pass_extra

