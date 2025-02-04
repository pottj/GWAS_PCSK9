#!/bin/bash

## Example SLURM script for BSU icelake jobs

## Section 1: SLURM Commands

## All SLURM commands must be placed at the start of the file
## Full documentation can be found here: https://slurm.schedmd.com/sbatch.html

## Enter a short name for the job, to be shown in SLURM output
#SBATCH -J UKB_PCSK9_regenie

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
module load ceuadmin/regenie/3.2.9

# If using the GPU cluster, replace the third line with the uncommented line:
# module load rhel8/default-amp

#! Insert additional module load commands after this line if needed:

## - - - - - - - - - - -

## Section 3: Run your application

## Step 1: create the necessary sample files for the --keep, --phenoFile and --covarFile

R CMD BATCH --vanilla ../scripts/01_Prep_01_checkUKB_samples.R ../scripts/01_Prep_01_checkUKB_samples.R.out
cp Rplots.pdf 01_Prep_01_Rplots.pdf
rm Rplots.pdf

## Step 2: preparing the genotype file
#
# Step 1 of regenie needs 1 genotype file
# This uses PLINK 1.9, not PLINK 2 - check if this is still necessary

rm -f /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/list_beds.txt
for chr in {2..22}; do echo "/rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c${chr}_b0_v2.bed /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb_snp_chr${chr}_v2.bim /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2_s488131.fam" >> /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/list_beds.txt; done

plink \
  --bed /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2.bed \
  --bim /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb_snp_chr1_v2.bim \
  --fam /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2_s488131.fam \
  --threads 20 \
  --keep /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_SampleList_all.txt \
  --merge-list /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/list_beds.txt \
  --make-bed --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/ukb_cal_allChrs
  
## Step 3: Exclusion files
#
# QC applied to the generated genotype file

module unload plink
module load plink/2.00-alpha

plink2 \
  --bfile /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/ukb_cal_allChrs \
  --keep /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_SampleList_White.txt \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --mind 0.1 \
  --threads 20 \
  --write-snplist --write-samples --no-id-header \
  --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/qc_pass
  
# 3. regenie step 1
#
# get LOCO

regenie \
  --step 1 \
  --bed /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/ukb_cal_allChrs \
  --extract /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/qc_pass.snplist \
  --keep /rds/user/jp2047/hpc-work/GWAS_PCSK9/UKB_genetics/qc_pass.id \
  --phenoFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_phenotypes_PCSK9.txt \
  --covarFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_covariates_noAncestry.txt \
  --threads 20 \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix /rds/user/jp2047/hpc-work/GWAS_PCSK9/regenie/tmpdir/regenie_tmp_preds \
  --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/regenie/ukb_step1_PCSK9
  
# 4. regenie step 2
#
# get association per chromosome, here test with chr 1

regenie \
  --step 2 \
  --bgen /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c1_b0_v3.bgen \
  --ref-first \
  --sample /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c1_b0_v3_s487160.sample \
  --phenoFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_phenotypes_PCSK9.txt \
  --covarFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/01_Prep_01_ukb_covariates_noAncestry.txt \
  --firth --approx --pThresh 0.01 \
  --pred /rds/user/jp2047/hpc-work/GWAS_PCSK9/regenie/ukb_step1_BT_pred.list \
  --bsize 400 \
  --threads 20 \
  --split \
  --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/regenie/ukb_step2_PCSK9_chr1
 

###############################################################
### You should not have to change anything below this line ####
###############################################################

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
if [ $SLURM_JOB_NUM_NODES -gt 1 ]; then
        echo "Running on nodes: $SLURM_JOB_NODELIST"
else
        echo "Running on node: `hostname`"
fi

echo "Current directory: `pwd`"
echo -e "\nNum tasks = $SLURM_NTASKS, Num nodes = $SLURM_JOB_NUM_NODES, OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD



