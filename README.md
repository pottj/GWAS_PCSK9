# GWAS of PCSK9 levels in the UK Biobank

last updated: 04/02/2025

## Aim 

Last year, I published my meta-GWAS of PCSK9, including 6 European cohorts (LIFE-Heart, LIFE-Adult, LURIC, TwinGene, KORA-F3, and GCKD). This was done double-stratified for both sex and lipid medication (using ATC Classification System code C10 for lipid modifying agents). 

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

1) Prepare PCSK9 data from the UKB 
2) Run GWAS using REGENIE
3) Run meta-GWAS for men, women, statin-free and statin-treated groups
4) Create data files with same structure as in meta-GWAS
5) GWAS replication check: test SNPs from [Pott et al., 2024](https://doi.org/10.1186/s13293-024-00602-6) if still genome-wide significantly associated
6) Run GCTA to identify independent signals at PCSK9 gene region
7) Interaction test: compare effect between 
    - pre- and post-menopausal women
    - men and post-menopausal women
    - men and pre-menopausal women
    - statin-free and statin-treated
8) Interaction replication check: test replicated signals for sex- and statin interaction
9) MR replication check: PCSK9 on LDL-C per subgroup
    - create LDL-C summary statistics with similar strata (UKB only)
    - test cis-MR only 
10) MVMR test for time-varying effect on LDL-C and/or CAD
    - using publicly available summary statistics
    - test cis instruments only
11) Combine meta-GWAS and UKB data 
    - women, statin-treated
    - women, statin-free
    - women
    - men, statin-treated
    - men, statin-free
    - men
    - sex-combined, statin-treated
    - sex-combined, statin-free
    - sex-combined, statin-combined

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
