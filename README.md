# GWAS of PCSK9 levels in the UK Biobank

last updated: 15/01/2025

## Aim 

In 2024, I published a meta-GWAS of PCSK9, including 6 European cohorts (LIFE-Heart, LIFE-Adult, LURIC, TwinGene, KORA-F3, and GCKD) ([Pott et al., 2024](https://doi.org/10.1186/s13293-024-00602-6)). This was done double-stratified for both sex and lipid medication (using ATC Classification System code C10 for lipid modifying agents). 

Now, I want to repeat this analysis in the UKB (replication approach). I want to replicate:

- the GWAS findings
- the observed SNP interactions (sex-interactions, statin-interaction)
- the observed causal effects of PCSK9 on LDL-C and their interactions 

In addition, I want to check if there are different effects between pre- and post-menopausal women, and if the causal effect on LDL-C or CAD is time-varying for women, using an MVMR approach. 

## Data sets

- UKB data from the MRC BSU (application 98032)
- PCSK9 summary statistics from my latest publication [Zenodo](https://zenodo.org/records/10600167)
- CAD summary statistics by sex [Agaram et al., 2022](https://hugeamp.org/dinspector.html?dataset=Aragam2022_CAD_EU) (EUR only) 
- LDL-C summary statistics by sex [Kanoni et al., 2022](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/sex_and_ancestry_specific_summary_stats/) (EUR only)

## Analysis plans

1) Prepare PCSK9 and LDL-C data from the UKB
    - Filter for pairwise kinship less than 10, consent still active (as of 18/08/2025), genetic sex == data base sex, no missing information for smoking or BMI
    - lipid lowering medication: ATC C10 coding
    - menopausal status: 
        - premenopausal: had menopause == no AND age<=50
        - postmenopausal: had menopause == yes AND age>=60
    - PCSK9 data set: anyone with PCSK9 measurement
    - LDL-C data set: anyone without PCSK9 measurement
    - 7 groups: 
        - men: sex == 1, with and without lipid lowering treatment
        - women: sex == 2, with and without lipid lowering treatment
        - postmenopausal: status == postmenopausal, with and without lipid lowering treatment
        - premenopausal: status == premenopausal, without lipid lowering treatment (only 200 treated)
2) Prepare UKB genetic data (PLINK 1.9)
    - create per biomarker one PLINK BED file set for REGENIE step 1
    - exclude SNPs with maf<0.01, mac<100, geno>=0.1 and HWE p-value<=1e-15
    - exclude samples for mind>= 0.1 (sample genotype missing rate)
3) Run REGENIE (v3.2.9) Step 1
    - done for each strata separately, as otherwise REGENIE would mean impute the missing phenotype data
    - adjusted for age, age squared, current smoking, log(BMI), genetic array (axiom or believe), genetic PCs (1-10)[, PPP batch number, PPP probe storage time on ice (time difference between blood draw and PPP measurement), in case of PCSK9]
4) Run REGENIE (v3.2.9) Step 2
    - done per chromosome (bgen files)
5) Create Summary Statistics file
    - same format as in previous work
6) Identify associated loci
7) Test if hits from previous work are replicated
8) Interaction test 
    - pre- and post-menopausal women
    - men and post-menopausal women
    - men and pre-menopausal women
    - statin-free and statin-treated
9) Mendelian Randomization of PCSK9 and LDL-C on CAD
10) Combine meta-GWAS and UKB data 

## Abbreviations

- ATC, Anatomical Therapeutic Chemical
- BSU, BioStatistics Unit
- CAD, Coronary Artery Disease
- GCTA, Genome-wide Complex Trait Analysis
- GWAS, Genome-Wide Association Study
- LDL-C, Low-Density Lipoprotein cholesterol
- MRC, Medical Research Council
- MVMR, MultiVariable Mendelian Randomization
- PCSK9, Proprotein Convertase Subtilisin/Kexin type 9
- SNP, Single Nucleotide Polymorphism
- UKB, UK Biobank
