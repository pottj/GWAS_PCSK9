## Step 3: run regenie step 1 per phenotype
#
# Whole genome regression model is fit at a subset of the total set of available genetic markers and get Leave One Chromosome Out (LOCO) predictions
mkdir -p /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/
mkdir -p /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/
mkdir -p /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/
mkdir -p /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/tmpdir_LDLC_${SLURM_ARRAY_TASK_ID}/
mkdir -p /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/tmpdir_PCSK9_${SLURM_ARRAY_TASK_ID}/

regenie \
  --step 1 \
  --bed /rds/user/jp2047/hpc-work/GWAS_PCSK9/03_UKB_genetics/ukb_PCSK9_set \
  --phenoFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_PCSK9_${SLURM_ARRAY_TASK_ID}.txt \
  --covarFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/02_regenie_input/01_Prep_ukb_covariates_EUR_PCSK9.txt \
  --catCovarList PPP_batch \
  --threads $SLURM_CPUS_PER_TASK \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/tmpdir_PCSK9_${SLURM_ARRAY_TASK_ID}/regenie_tmp_preds \
  --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/pheno_${SLURM_ARRAY_TASK_ID}

regenie \
  --step 1 \
  --bed /rds/user/jp2047/hpc-work/GWAS_PCSK9/03_UKB_genetics/ukb_LDLC_set \
  --phenoFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_LDLC_${SLURM_ARRAY_TASK_ID}.txt \
  --covarFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/02_regenie_input/01_Prep_ukb_covariates_EUR_LDLC.txt \
  --threads $SLURM_CPUS_PER_TASK \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/tmpdir_LDLC_${SLURM_ARRAY_TASK_ID}/regenie_tmp_preds \
  --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/pheno_${SLURM_ARRAY_TASK_ID}

cat /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/pheno_1_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/pheno_2_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/pheno_3_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/pheno_4_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/pheno_5_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/pheno_6_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/pheno_7_pred.list > /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/pheno_pred.list
  
cat /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/pheno_1_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/pheno_2_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/pheno_3_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/pheno_4_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/pheno_5_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/pheno_6_pred.list /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/pheno_7_pred.list > /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/pheno_pred.list

rm -r /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/tmpdir_LDLC_${SLURM_ARRAY_TASK_ID}/
rm -r /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/tmpdir_PCSK9_${SLURM_ARRAY_TASK_ID}/
  
  