## Step 2: preparing the genotype file
#
# Step 1 of regenie needs 1 genotype file (not just imputed genedosages)
mkdir -p /rds/user/jp2047/hpc-work/GWAS_PCSK9/03_UKB_genetics/
rm -f /rds/user/jp2047/hpc-work/GWAS_PCSK9/03_UKB_genetics/list_beds_2_to_22.txt

echo Running on $SLURM_CPUS_PER_TASK cores

for chr in {2..22}; do echo "/rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c${chr}_b0_v2.bed /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb_snp_chr${chr}_v2.bim /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2_s488131.fam" >> /rds/user/jp2047/hpc-work/GWAS_PCSK9/03_UKB_genetics/list_beds_2_to_22.txt; done

# Create plink files for PCSK9 run
plink \
  --bed /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2.bed \
  --bim /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb_snp_chr1_v2.bim \
  --fam /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2_s488131.fam \
  --threads $SLURM_CPUS_PER_TASK \
  --keep /rds/user/jp2047/hpc-work/GWAS_PCSK9/02_regenie_input/01_Prep_ukb_SampleList_EUR_PCSK9.txt \
  --merge-list /rds/user/jp2047/hpc-work/GWAS_PCSK9/03_UKB_genetics/list_beds_2_to_22.txt \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --mind 0.1 \
  --make-bed --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/03_UKB_genetics/ukb_PCSK9_set

# Create plink files for PCSK9 run
plink \
  --bed /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2.bed \
  --bim /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb_snp_chr1_v2.bim \
  --fam /rds/user/jp2047/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c1_b0_v2_s488131.fam \
  --threads $SLURM_CPUS_PER_TASK \
  --keep /rds/user/jp2047/hpc-work/GWAS_PCSK9/02_regenie_input/01_Prep_ukb_SampleList_EUR_LDLC.txt \
  --merge-list /rds/user/jp2047/hpc-work/GWAS_PCSK9/03_UKB_genetics/list_beds_2_to_22.txt \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --mind 0.1 \
  --make-bed --out /rds/user/jp2047/hpc-work/GWAS_PCSK9/03_UKB_genetics/ukb_LDLC_set


# Check if exit code != 0 for last command (plink)
retVal =$?
if [ $retVal -ne 0 ]; then
    echo "Error"
    exit $retVal 
fi 


