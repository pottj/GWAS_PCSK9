## Step 4: run regenie step 2 per phenotype
#
mkdir -p /rds/user/jp2047/hpc-work/GWAS_PCSK9/05_regenie_step2/
mkdir -p /rds/user/jp2047/hpc-work/GWAS_PCSK9/05_regenie_step2/PCSK9/
mkdir -p /rds/user/jp2047/hpc-work/GWAS_PCSK9/05_regenie_step2/LDLC/

regenie \
  --step 2 \
  --bgen /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c${SLURM_ARRAY_TASK_ID}_b0_v3.bgen \
  --ref-first \
  --sample /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c${SLURM_ARRAY_TASK_ID}_b0_v3_s487160.sample \
  --phenoFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step2_PCSK9.txt \
  --covarFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/02_regenie_input/01_Prep_ukb_covariates_EUR_PCSK9.txt \
  --catCovarList PPP_batch \
  --firth --approx --pThresh 0.01 \
  --pred /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/PCSK9/pheno_pred.list \
  --bsize 400 \
  --threads $SLURM_CPUS_PER_TASK \
  --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/05_regenie_step2/PCSK9/ukb_chr${SLURM_ARRAY_TASK_ID}

regenie \
  --step 2 \
  --bgen /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c${SLURM_ARRAY_TASK_ID}_b0_v3.bgen \
  --ref-first \
  --sample /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c${SLURM_ARRAY_TASK_ID}_b0_v3_s487160.sample \
  --phenoFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step2_LDLC.txt \
  --covarFile /rds/user/jp2047/hpc-work/GWAS_PCSK9/02_regenie_input/01_Prep_ukb_covariates_EUR_LDLC.txt \
  --firth --approx --pThresh 0.01 \
  --pred /rds/user/jp2047/hpc-work/GWAS_PCSK9/04_regenie_step1/LDLC/pheno_pred.list \
  --bsize 400 \
  --threads $SLURM_CPUS_PER_TASK \
  --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/05_regenie_step2/LDLC/ukb_LDLC_chr${SLURM_ARRAY_TASK_ID}


