#' ---
#' title: "Check UKB data"
#' subtitle: "GWAS PCSK9 (Olink)"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")

#' # Load UKB Data ####
#' ***
{
  #' I want the proteomics data from the BSU UKB application and all relevant covariables: 
  #' 
  #' Field ID   | Description
  #' 
  #' eid        | Encoded anonymised participant ID
  #'    31      | Sex (0 = female, 1 = male)
  #'    53      | date of measurement
  #'    54      | UKB center  
  #'  2724      | Had menopause (0 = no, 1 = yes)
  #'  3700      | Time since last menstrual period
  #'  3710      | Length of menstrual cycle
  #'  3720      | Menstruating today (0 = no, 1 = yes)
  #' 20003      | Treatment/medication code (48 columns)
  #' 20116      | Smoking status (0 = never, 1 = previous, 2 = current)
  #' 21000      | Ethnic background 
  #' 21001      | Body mass index (BMI)
  #' 21022      | Age at recruitment
  #' 22000      | Genotype measurement batch
  #' 22001      | Genetic sex
  #' 22006      | Genetic ethnic grouping
  #' 22009      | Genetic principal components (40 columns)
  #' 22021      | Genetic kinship to other participants
  #' 30780      | LDL direct
  #' 30800      | Estradiol
  #' 30850      | Testosterone
  #' 
  myTab_20 = fread(UKB_phenotypes, header=TRUE, sep="\t",nrows = 20)
  myAnnot = data.table(colNm = names(myTab_20))
  myAnnot[,colNR := 1:18506]
  
  covars = c("f.eid", "f.31.0.0","f.53.0.0","f.54.0.0", 
             "f.2724.0.0","f.3700.0.0","f.3710.0.0","f.3720.0.0", 
             paste0("f.20003.0.",c(0:47)),"f.20116.0.0","f.21000.0.0","f.21001.0.0","f.21022.0.0",
             "f.22000.0.0","f.22001.0.0","f.22006.0.0",paste0("f.22009.0.",1:10), "f.22021.0.0", 
             "f.30780.0.0","f.30800.0.0","f.30850.0.0")
  table(is.element(covars,myAnnot$colNm))
  
  myAnnot = myAnnot[colNm %in% c(covars),]
  
  x = myAnnot[,colNR]
  myTab = fread(UKB_phenotypes, header=TRUE, sep="\t",select = x)
  names(myTab)
  names(myTab) = c("ID","sex","date","center",
                   "menopause","lastMenst","cycleLength","todayMenst",
                   paste("meds",1:48,sep="_"),"smoking","ancestry","BMI","age",
                   "genetic_batch","genetic_sex","genetic_ancestry",paste("PC",1:10,sep="_"),"kinship", 
                   "LDLC","E2","TT")
  
  #' Save data
  save(myTab, file = paste0(data_QC,"/01_UKB/01_Prep_a_unfiltered.RData"))
  #load(paste0(data_QC,"/01_UKB/01_Prep_a_unfiltered.RData"))
}

#' # Data QC #### 
#' ***
#' ## Hard Filtering
#' 
#' - Pairwise kinship less than 10
#' - Consent still active
#' - Genetic sex = data base sex
#' - Missing information for smoking or BMI
#' 
{
  #' Filter for people without any kinship
  myTab[,table(kinship)]
  myTab[,table(is.na(kinship))]
  myTab = myTab[!is.na(kinship),]
  myTab = myTab[kinship != 10,]
  
  #' Filter for consent
  ToExclude = fread(gsub("ukb672224.tab","withdraw98032_20250818.csv",UKB_phenotypes))
  table(is.element(myTab$ID,ToExclude$V1))
  myTab = myTab[!is.element(ID,ToExclude$V1),]
  
  #' Filter for same sex information
  myTab[,table(sex,genetic_sex)]
  myTab = myTab[sex == genetic_sex,]
  myTab[sex==0,sex := 2]
  myTab[,genetic_sex := NULL]
  
  #' Check smoking: I want active smoking coded as 1, and never and ex-smoking coded as 0
  myTab[,table(is.na(smoking))]
  myTab = myTab[!is.na(smoking)]
  myTab[,table(smoking)]
  myTab = myTab[smoking != -3,]
  myTab[smoking == 1,smoking := 0]
  myTab[smoking == 2,smoking := 1]
  
  #' Check BMI: how is it distributed
  myTab[,hist(BMI)]
  myTab[,table(BMI>50,BMI<15)]
  myTab = myTab[!is.na(BMI) & BMI>=15 & BMI<=50,]
  
  #' Check dimension after filtering and save
  dim(myTab)
  save(myTab, file = paste0(data_QC,"/01_UKB/01_Prep_b_filtered.RData"))
  #load(paste0(data_QC,"/01_UKB/01_Prep_b_filtered.RData"))
  
}

#' ## New Variables 
#' 
#' - lipid lowering medication: ATC C10, binary
#' - assessment center: city names
#' - genotyping batch: switch to array type (Believe, Axiom)
#' - genetic ancestry/ethnicity: merge groups
#' - menopause: create groups of pre- and post-menopausal women
#' 
{
  #' **Check medication**: I want information on lipid lowering medication (ATC starting with C10, Lipid modifying agents). I will use the supplemental data 1 file of Wu et al. (Wu, Y., Byrne, E.M., Zheng, Z. et al. Genome-wide association study of medication-use and associated disease in the UK Biobank. Nat Commun 10, 1891 (2019). https://doi.org/10.1038/s41467-019-09572-5)
  #' 
  MedCoding = data.table(read_xlsx(UKB_MedicationCoding,sheet=1))
  lipidMedication = MedCoding[grepl("C10",ATC_Code),Coding]
  
  myMeds = names(myTab)[grep("meds",names(myTab))]
  myTab[,lipidMeds := 0]
  for(i in 1:length(myMeds)){
    #i=1
    myTab[get(myMeds[i]) %in% lipidMedication,lipidMeds := 1]
  }
  myTab[,get("myMeds"):=NULL]
  table(myTab$lipidMeds)
 
  #' **Check assessment center**: I want to switch from the numbers to characters, so I do not by accident use them as normal variable but as factors. 
  myTab[,table(center)]
  CenterCoding = fread(UKB_CenterCoding)
  matched = match(myTab$center,CenterCoding$coding)
  table(is.na(matched))
  myTab[,centerName := CenterCoding[matched,meaning]]
  myTab[,table(centerName)]
  
  #' **Check genotyping batch**: There are potential batch effects due to the different array used for genotyping. I to not need the exact batch number, but the array type. 
  myTab[,table(genetic_batch)]
  myTab[,genetic_array := "Axiom"]
  myTab[genetic_batch<0,genetic_array := "BiLEVE"]
  myTab[,table(genetic_array, sex)]
  myTab[,genetic_batch := NULL]
  
  #' **Check genetic ancestry**: for the GWAS, I will only use samples who self-identified as 'White British' according to Field 21000 and have very similar genetic ancestry based on a principal components analysis of the genotypes (--> genetic ancestry). For plotting & comparison reason, I keep the ethnic background.  
  myTab[,table(is.na(genetic_ancestry),ancestry)]
  
  #' **Check ethnicity**:
  myTab[,table(ancestry)]
  myTab[,ancestryGroups := substr(ancestry,1,1)]
  myTab[ancestryGroups == "-",ancestryGroups := NA]
  myTab[,table(ancestryGroups,genetic_ancestry)]
  myTab[,table(ancestry,ancestryGroups)]
  
  #' **Check sex groups**: 
  #' - pre-menopausal: all women who stated they are still menstruating and who are younger than 51
  #' - post-menopausal: all women who stated they are no longer menstruating and who are older than 60
  myTab[sex==2 & (menopause==1 & age>=60),group := "postmenopausal"]
  myTab[sex==2 & (menopause==0 & age<=50),group := "premenopausal"]
  myTab[,table(group)]
  myTab[,table(is.na(group),menopause)]
  
  #' Check dimension after filtering and save
  dim(myTab)
  save(myTab, file = paste0(data_QC,"/01_UKB/01_Prep_c_filtered.RData"))
  #load(paste0(data_QC,"/01_UKB/01_Prep_c_filtered.RData"))
  
}

#' ## Proteomics data
#' 
#' - PCSK9 measurement
#' - PPP sample information: ukb676343.csv (data field 30900 - 30903)
#' - PPP batch: olink_batch_number.dat (see https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=1839, Resources), 
#' - PPP time between blood sampling and measurement: olink_processing_start_date.dat (see https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=1839, Resources)

{
  olink_data = fread(UKB_proteomics_data)
  olink_coding = fread(UKB_proteomics_coding)
  olink_coding[,gene:= gsub(";.*","",meaning)]
  olink_coding[,description:= gsub(".*;","",meaning)]
  olink_data = olink_data[protein_id == olink_coding[gene=="PCSK9",coding] & ins_index==0,]
  
  matched = match(myTab$ID,olink_data$eid)
  myTab[,PCSK9 := olink_data[matched,result]]
  
  olink_samples = fread(paste0(UKB_proteomics,"/ukb678244.csv"))
  olink_samples = olink_samples[!is.na(`30900-0.0`),]
  
  matched = match(myTab$ID,olink_samples$eid)
  myTab[,PPP_NRprot := olink_samples[matched,`30900-0.0`],]
  myTab[,PPP_plate := olink_samples[matched,`30901-0.0`],]
  myTab[,PPP_well := olink_samples[matched,`30902-0.0`],]
  
  olink_batch = fread(paste0(UKB_proteomics,"/olink_batch_number.dat"))
  olink_date = fread(paste0(UKB_proteomics,"/olink_processing_start_date.dat"))
  matched = match(olink_date$PlateID,olink_batch$PlateID)
  olink_date[,batch := olink_batch[matched,Batch]]
  
  matched = match(myTab$PPP_plate,olink_date$PlateID)
  table(is.na(matched))
  myTab[,PPP_batch := olink_date[matched,batch]]
  myTab[,PPP_Startdate := olink_date[matched,Processing_StartDate]]
  myTab[,PPP_timeDif := (PPP_Startdate - date)/365.25]
  myTab[,PPP_timeDif := as.numeric(PPP_timeDif)]
  
  myTab[,PPP_plate := gsub("890000000","plate_",PPP_plate)]
  myTab[,PPP_batch := paste0("batch_",PPP_batch)]
  myTab[,table(PPP_batch)]
  
  # everyone without PCSK9 measurement should be NA
  myTab[is.na(PCSK9),PPP_NRprot := NA]
  myTab[is.na(PCSK9),PPP_plate := NA]
  myTab[is.na(PCSK9),PPP_well := NA]
  myTab[is.na(PCSK9),PPP_batch := NA]
  myTab[is.na(PCSK9),PPP_Startdate := NA]
  myTab[is.na(PCSK9),PPP_timeDif := NA]

  #' Check dimension after filtering and save
  dim(myTab)
  save(myTab, file = paste0(data_QC,"/01_UKB/01_Prep_d_filtered.RData"))
  #load(paste0(data_QC,"/01_UKB/01_Prep_c_filtered.RData"))
  
}

#' # Check non-genetic regression ####
#' ***
myTab[,age2 := age^2]

summary(lm(PCSK9 ~ sex + smoking + age + age2 + lipidMeds + log(BMI), data = myTab))

summary(lm(PCSK9 ~ smoking + age + age2 + lipidMeds + log(BMI), data = myTab, subset = sex==1))
summary(lm(PCSK9 ~ smoking + age + age2 + lipidMeds + log(BMI), data = myTab, subset = sex==2))

summary(lm(PCSK9 ~ smoking + age + age2 + lipidMeds + log(BMI), data = myTab, subset = group=="premenopausal"))
summary(lm(PCSK9 ~ smoking + age + age2 + lipidMeds + log(BMI), data = myTab, subset = group=="postmenopausal"))

summary(lm(PCSK9 ~ smoking + age + age2 + lipidMeds + log(BMI) + TT*sex + as.factor(ancestryGroups) + PPP_timeDif, data = myTab))

#' **Conclusion**
#' 
#' I will keep using the same model as in my previous meta-GWAS, but will run: 
#' 
#' - men per treatment groups
#' - women per treatment groups
#' - premenopausal women treatment free only (not enough treated women)
#' - postmenopausal women per treatment groups
#' 
#' Check dimension after filtering and save
dim(myTab)
save(myTab, file = paste0(data_QC,"/01_UKB/01_Prep_e_final.RData"))
#load(paste0(data_QC,"/01_UKB/01_Prep_e_final.RData"))

#' # Plots ####
#' ***
#' - boxplot per group
#' - boxplot per ancestry
#' - boxplot per lipid lowering medication and group
#' - scatterplot per BMI and group
#' 
plotData = copy(myTab)
plotData[,sex2 := "men"]
plotData[sex==2,sex2 := "women"]
plotData[,lipidMeds2 := "yes"]
plotData[lipidMeds==0,lipidMeds2 := "no"]
plotData[,strata := paste0(sex2," - ",lipidMeds2)]
plotData[,strata2 := paste0(group," - ",lipidMeds2)]
plotData[lipidMeds==0,strata3 := paste0(sex2," - no medication")]
plotData[lipidMeds==1,strata3 := paste0(sex2," - medication")]
plotData[lipidMeds==0,strata4 := paste0(group," - no medication")]
plotData[lipidMeds==1,strata4 := paste0(group," - medication")]

plotData[,smoking2 := "1 - current"]
plotData[smoking==0,smoking2 := "0 - never & previous"]
plotData[ancestry< -2,ancestry2 := "7 - prefer not to answer"]
plotData[ancestry< 0 & ancestry >-2,ancestry2 := "8 - do not know"]
plotData[(ancestry< 1005 & ancestry > 10) | ancestry == 1,ancestry2 := "1 - White"]
plotData[(ancestry< 2005 & ancestry > 1005) | ancestry == 2,ancestry2 := "2 - Mixed"]
plotData[(ancestry< 3005 & ancestry > 2005) | ancestry == 3,ancestry2 := "3 - Asian (British)"]
plotData[(ancestry< 4005 & ancestry > 3005) | ancestry == 4,ancestry2 := "4 - Black (British)"]
plotData[ancestry == 5,ancestry2 := "5 - Chinese"]
plotData[ancestry == 6,ancestry2 := "6 - Other"]
plotData[,ancestry3 := ancestry]
plotData[ancestry %in% c(1,2,3,4,5,6),ancestry3 := ancestry * 1000]
plotData[,min(age)]
plotData[,max(age)]
plotData[age<=45, age3 := "40-45"]
plotData[age<=50 & age>45, age3 := "46-50"]
plotData[age<=55 & age>50, age3 := "51-55"]
plotData[age<=60 & age>55, age3 := "56-60"]
plotData[age<=65 & age>60, age3 := "61-65"]
plotData[age<=70 & age>65, age3 := "66-70"]

p1.1 = ggplot(data = plotData[!is.na(PCSK9)], aes(x = lipidMeds2, y = PCSK9, fill=strata)) + 
  facet_wrap(~sex2) + 
  geom_boxplot() +
  labs(x="lipid lowering medication",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels by sex and lipid-lowering medication")) +
  scale_fill_manual(values = c("#0F9ED5","#CCDFEF","#E97132","#F7D5CD")) +
  theme_bw() +
  theme(legend.position = "none")
p1.1
ggsave(plot = p1.1, filename = paste0("../results/_figures/01_PrepData/Boxplot_LipLowMed_Sex.png"), 
       height = 7, width = 14)

p1.2 = ggplot(data = plotData[!is.na(PCSK9) & !is.na(group)], aes(x = lipidMeds2, y = PCSK9, fill=strata2)) + 
  facet_wrap(~group) + 
  geom_boxplot() +
  labs(x="lipid lowering medication",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels by menopausal status and lipid-lowering medication")) +
  scale_fill_manual(values = c("#4EA72E","#D0E1CD","#A02B93","#DFCDDC")) +
  theme_bw() +
  theme(legend.position = "none")
p1.2
ggsave(plot = p1.2, filename = paste0("../results/_figures/01_PrepData/Boxplot_LipLowMed_Menopause.png"), 
       height = 7, width = 14)

p2.1 = ggplot(plotData[!is.na(PCSK9),],aes(x=PCSK9,y=log(BMI))) + 
  facet_wrap(~ strata3,scales = "fixed") + 
  #geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
  geom_point() + 
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "pearson") + 
  xlab("PCSK9 levels") + 
  ylab("log(BMI)")+
  theme_bw()
p2.1
ggsave(plot = p2.1, filename = paste0("../results/_figures/01_PrepData/Scatterplot_BMI_Strata.png"), 
       height = 7, width = 14)

p2.2 = ggplot(plotData[!is.na(PCSK9) & !is.na(group),],aes(x=PCSK9,y=log(BMI))) + 
  facet_wrap(~ strata4,scales = "fixed") + 
  #geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
  geom_point() + 
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "pearson") + 
  xlab("PCSK9 levels") + 
  ylab("log(BMI)")+
  theme_bw()
p2.2
ggsave(plot = p2.2, filename = paste0("../results/_figures/01_PrepData/Scatterplot_BMI_Menopause.png"), 
       height = 7, width = 14)

p3.1 = ggplot(plotData[!is.na(PCSK9) & !is.na(LDLC),],aes(x=PCSK9,y=LDLC)) + 
  facet_wrap(~ strata3,scales = "fixed") + 
  #geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
  geom_point() + 
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "pearson") + 
  xlab("PCSK9 levels") + 
  ylab("LDL-C levels")+
  theme_bw()
p3.1
ggsave(plot = p3.1, filename = paste0("../results/_figures/01_PrepData/Scatterplot_LDLC_Strata.png"), 
       height = 7, width = 14)

p3.2 = ggplot(plotData[!is.na(PCSK9) & !is.na(group) & !is.na(LDLC),],aes(x=PCSK9,y=LDLC)) + 
  facet_wrap(~ strata4,scales = "fixed") + 
  #geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
  geom_point() + 
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "pearson") + 
  xlab("PCSK9 levels") + 
  ylab("LDL-C levels")+
  theme_bw()
p3.2
ggsave(plot = p3.2, filename = paste0("../results/_figures/01_PrepData/Scatterplot_LDLC_Menopause.png"), 
       height = 7, width = 14)

p4 = ggplot(data = plotData[!is.na(PCSK9)], aes(x = as.factor(ancestry2), y = PCSK9, fill = as.factor(ancestry3))) + 
  geom_boxplot() +
  labs(x="Ethnic background",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels by ethnic background")) +
  theme_bw() +
  theme(legend.position = "none")
p4
ggsave(plot = p4, filename = paste0("../results/_figures/01_PrepData/Boxplot_Ancestry.png"), 
       height = 7, width = 14)

p5.1 = ggplot(data = plotData[!is.na(PCSK9)], aes(x = strata3, y = PCSK9, fill = strata)) + 
  facet_wrap(~smoking2) + 
  geom_boxplot() +
  labs(x="current smoking",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels by current smoking and strata")) +
  scale_fill_manual(values = c("#0F9ED5","#CCDFEF","#E97132","#F7D5CD")) +
  theme_bw() +
  theme(legend.position = "none")
p5.1
ggsave(plot = p5.1, filename = paste0("../results/_figures/01_PrepData/Boxplot_Smoking_Strata.png"), 
       height = 7, width = 14)

p5.2 = ggplot(data = plotData[!is.na(PCSK9) & !is.na(group)], aes(x = strata4, y = PCSK9, fill = strata2)) + 
  facet_wrap(~smoking2) + 
  geom_boxplot() +
  labs(x="current smoking",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels by current smoking and menopausal status")) +
  scale_fill_manual(values = c("#4EA72E","#D0E1CD","#A02B93","#DFCDDC")) +
  theme_bw() +
  theme(legend.position = "none")
p5.2
ggsave(plot = p5.2, filename = paste0("../results/_figures/01_PrepData/Boxplot_Smoking_Menopause.png"), 
       height = 7, width = 14)

p6.1 = ggplot(data = plotData[sex==1 & !is.na(PCSK9)], aes(x = age3, y = PCSK9, fill = lipidMeds2)) + 
  #facet_wrap(~lipidMeds2) + 
  geom_boxplot() +
  labs(x="age (in years)",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels and age by lipid lowering medication in men")) +
  scale_fill_manual(values = c("#0F9ED5","#CCDFEF"),
                    labels = c("no", "yes"))+
  theme_bw() +
  theme(legend.position = "none")
p6.1
ggsave(plot = p6.1, filename = paste0("../results/_figures/01_PrepData/Boxplot_Age_Men.png"), 
       height = 7, width = 14)

p6.2 = ggplot(data = plotData[sex==2 & !is.na(PCSK9)], aes(x = age3, y = PCSK9, fill = lipidMeds2)) + 
  #facet_wrap(~lipidMeds2) + 
  geom_boxplot() +
  labs(x="age (in years)",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels and age by lipid lowering medication in women")) +
  scale_fill_manual(values = c("#E97132","#F7D5CD"),
                    labels = c("no", "yes"))+
  theme_bw() +
  theme(legend.position = "none")
p6.2
ggsave(plot = p6.2, filename = paste0("../results/_figures/01_PrepData/Boxplot_Age_Women.png"), 
       height = 7, width = 14)

#' # Create input for regenie ####
#' ***
#' ## Filter for EUR ancestry
#' 
myTab = myTab[genetic_ancestry == 1,]

#' Check dimension after filtering and save
dim(myTab)
save(myTab, file = paste0(data_QC,"/01_UKB/01_Prep_e_final_EUR.RData"))

#' ## Sample inclusion 
#' 
#' **Sample inclusion/exclusion file format**: No header. Each line starts with individual FID IID. Space/tab separated.
#' 
#' Needed for the generation of the bed/bim/fam files for regenie step 1
#' 
myTab1 = copy(myTab)
myTab1[,IID := ID]
setnames(myTab1,"ID","FID")
filt_PCSK9 = !is.na(myTab1$PCSK9)
myTab1 = myTab1[,c(1,41)]
write.table(myTab1[filt_PCSK9,],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_SampleList_EUR_PCSK9.txt"), 
            col.names = F, row.names = F, quote = F)
write.table(myTab1[!filt_PCSK9,],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_SampleList_EUR_LDLC.txt"), 
            col.names = F, row.names = F, quote = F)

#' ## Covariate file 
#' 
#' **Covariate file format**: Line 1 : Header with FID, IID and C covariate names. Followed by lines of C+2 values. Space/tab separated. Each line contains individual FID and IID followed by C covariate values.
#' 
#' Same file for regenie step 1 and step 2.
#' 
#' This includes: age, smoking, BMI, gArray, gPC1 - 10
#' 
myNames = c("ID", "age", "age2", "smoking", "BMI", "genetic_array", paste0("PC_",1:10), "PPP_batch", "PPP_timeDif" ) 
myTab2 = copy(myTab)
colsOut<-setdiff(colnames(myTab2),myNames)
myTab2[,get("colsOut"):=NULL]
setcolorder(myTab2,myNames)
dim(myTab2)
myTab2 = myTab2[!is.na(PPP_batch)]
myTab2[,BMI := log(BMI)]
myTab2[,genetic_array := gsub("Axiom","1",genetic_array)]
myTab2[,genetic_array := gsub("BiLEVE","0",genetic_array)]
myTab2 = cbind(myTab2$ID,myTab2)
names(myTab2)[1:2] = c("FID","IID")

write.table(myTab2,
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_covariates_EUR_PCSK9.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")

#' Now for LDLC
myNames = c("ID", "age", "age2", "smoking", "BMI", "genetic_array", paste0("PC_",1:10),"PPP_batch") 
myTab3 = copy(myTab)
colsOut<-setdiff(colnames(myTab3),myNames)
myTab3[,get("colsOut"):=NULL]
setcolorder(myTab3,myNames)
dim(myTab3)
myTab3 = myTab3[is.na(PPP_batch)]
myTab3[,PPP_batch := NULL]
myTab3[,BMI := log(BMI)]
myTab3[,genetic_array := gsub("Axiom","1",genetic_array)]
myTab3[,genetic_array := gsub("BiLEVE","0",genetic_array)]
myTab3 = cbind(myTab3$ID,myTab3)
names(myTab3)[1:2] = c("FID","IID")

write.table(myTab3,
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_covariates_EUR_LDLC.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")

#' ## Phenotype file
#' 
#' **Phenotype file format**: Line 1 : Header with FID, IID and P phenotypes names. Followed by lines of P+2 values. Space/tab separated. Each line contains individual FID and IID followed by P phenotype values (for binary traits, must be coded as 0=control, 1=case, NA=missing unless using --1).
#' 
#' 7 files with one phenotype each for regenie step 1; and 1 file with all six phenotypes for regenie step 2. Reason: step 1 will mean-impute missing phenotypes, but I don't want that for my subgroups as they are independent (no overlaps).  
#' 
myTab3 = copy(myTab)
myTab3[sex==1 & lipidMeds==0,PCSK9_men_free := PCSK9]
myTab3[sex==1 & lipidMeds==1,PCSK9_men_treated := PCSK9]
myTab3[sex==2 & lipidMeds==0,PCSK9_women_free := PCSK9]
myTab3[sex==2 & lipidMeds==1,PCSK9_women_treated := PCSK9]
myTab3[group=="premenopausal" & lipidMeds==0,PCSK9_pre_free := PCSK9]
myTab3[group=="postmenopausal" & lipidMeds==0,PCSK9_post_free := PCSK9]
myTab3[group=="postmenopausal" & lipidMeds==1,PCSK9_post_treated := PCSK9]

myTab3[!is.na(PCSK9),LDLC := NA]
myTab3[sex==1 & lipidMeds==0,LDLC_men_free := LDLC]
myTab3[sex==1 & lipidMeds==1,LDLC_men_treated := LDLC]
myTab3[sex==2 & lipidMeds==0,LDLC_women_free := LDLC]
myTab3[sex==2 & lipidMeds==1,LDLC_women_treated := LDLC]
myTab3[group=="premenopausal" & lipidMeds==0,LDLC_pre_free := LDLC]
myTab3[group=="postmenopausal" & lipidMeds==0,LDLC_post_free := LDLC]
myTab3[group=="postmenopausal" & lipidMeds==1,LDLC_post_treated := LDLC]

myTab4 = cbind(myTab3$ID,myTab3$ID,myTab3[,c(41:47)])
names(myTab4)[1:2] = c("FID","IID")
filt_PCSK9 = !is.na(myTab3$PCSK9)
myTab4 = myTab4[filt_PCSK9,]

write.table(myTab4,
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step2_PCSK9.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab4[,c(1,2,3)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_PCSK9_1.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab4[,c(1,2,4)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_PCSK9_2.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab4[,c(1,2,5)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_PCSK9_3.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab4[,c(1,2,6)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_PCSK9_4.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab4[,c(1,2,7)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_PCSK9_5.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab4[,c(1,2,8)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_PCSK9_6.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab4[,c(1,2,9)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_PCSK9_7.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")

myTab5 = cbind(myTab3$ID,myTab3$ID,myTab3[,c(48:54)])
names(myTab5)[1:2] = c("FID","IID")
filt_PCSK9 = !is.na(myTab3$PCSK9)
myTab5 = myTab5[!filt_PCSK9,]

write.table(myTab5,
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step2_LDLC.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab5[,c(1,2,3)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_LDLC_1.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab5[,c(1,2,4)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_LDLC_2.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab5[,c(1,2,5)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_LDLC_3.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab5[,c(1,2,6)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_LDLC_4.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab5[,c(1,2,7)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_LDLC_5.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab5[,c(1,2,8)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_LDLC_6.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab5[,c(1,2,9)],
            file = paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step1_LDLC_7.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

