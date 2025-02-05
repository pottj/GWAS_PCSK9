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

#' # UKB Data ####
#' ***
#' I want the proteomics data from the BSU UKB application and all relevant covariables: 
#' 
#' ## Load all covariables ####
#' 
#' Field ID   | Description
#' 
#' eid        | Encoded anonymised participant ID
#'    31      | Sex (0 = female, 1 = male)
#'  2724      | Had menopause (0 = no, 1 = yes)
#'  3700      | Time since last menstrual period
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
#' 
myTab_20 = fread(UKB_phenotypes, header=TRUE, sep="\t",nrows = 20)
myAnnot = data.table(colNm = names(myTab_20))
myAnnot[,colNR := 1:18506]

covars = c("f.eid", "f.31.0.0", "f.2724.0.0","f.3700.0.0","f.3720.0.0", 
           paste0("f.20003.0.",c(0:47)),"f.20116.0.0",
           "f.21000.0.0","f.21001.0.0","f.21022.0.0","f.22000.0.0","f.22001.0.0","f.22006.0.0",
           paste0("f.22009.0.",1:10), "f.22021.0.0")
table(is.element(covars,myAnnot$colNm))

myAnnot = myAnnot[colNm %in% c(covars),]

x = myAnnot[,colNR]
myTab = fread(UKB_phenotypes, header=TRUE, sep="\t",select = x)
names(myTab)
names(myTab) = c("ID","sex","menopause","lastMenst","todayMenst",
                 paste("medication_init",1:48,sep="_"),"smoking",
                 "ancestry","BMI","age","gBatch","gSex","gAncestry",
                 paste("gPC",1:10,sep="_"),"gKinship")

#' Check medication: I will exclude everyone taking lipid-lowering medication (ATC starting with C10)
codingTable = data.table(read_xlsx(MedicationCoding,sheet=1))
medicationC10 = codingTable[grepl("C10",ATC_Code),Coding]

myMeds = names(myTab)[grep("medication_init",names(myTab))]
myTab[,lipidMeds := 0]
for(i in 1:length(myMeds)){
  #i=1
  myTab[get(myMeds[i]) %in% medicationC10,lipidMeds := 1]
}
myTab[,get("myMeds"):=NULL]
myTab[,table(lipidMeds,sex)]

#' Check for consent: There might be people who have withdrawn their consent by now. Last information on consent: 17/12/2024
ToExclude = fread(gsub("ukb672224.tab","withdraw98032_20241217.csv",UKB_phenotypes))
table(is.element(myTab$ID,ToExclude$V1))
myTab = myTab[!is.element(ID,ToExclude$V1),]

#' Check genotyping batch: There are potential batch effects due to the different array used for genotyping. I to not need the exact batch number, but the array type. 
myTab[,table(gBatch)]
myTab[,gArray := "Axiom"]
myTab[gBatch<0,gArray := "BiLEVE"]
myTab[,table(gArray, sex)]

#' Save as temporary file 
save(myTab, file = paste0(data_QC,"/01_Prep_01_ukb_PCSK9.RData"))
# load(paste0(data_QC,"/01_Prep_01_ukb_PCSK9.RData"))

#' ## Proteomics data ####
#' 
olink_samples = fread(UKB_proteomics_samples)
olink_samples = olink_samples[!is.na(`30900-0.0`),]
myTab = myTab[ID %in% olink_samples$eid,]

olink_data = fread(UKB_proteomics_data)
olink_coding = fread(UKB_proteomics_coding)
olink_coding[,gene:= gsub(";.*","",meaning)]
olink_coding[,description:= gsub(".*;","",meaning)]
olink_data = olink_data[protein_id == olink_coding[gene=="PCSK9",coding] & ins_index==0,]

matched = match(myTab$ID,olink_data$eid)
myTab[,PCSK9 := olink_data[matched,result]]
myTab = myTab[!is.na(PCSK9),]

#' Save as temporary file 
save(myTab, file = paste0(data_QC,"/01_Prep_01_ukb_PCSK9.RData"))
# load(paste0(data_QC,"/01_Prep_01_ukb_PCSK9.RData"))

#' # Checks ####
#' ***
#' ## Define Groups
#' I want the genetic sex matching the menopause information. So I will exclude all samples with mismatching genetic and data base sex. 
#' 
#' I also want only pre- and not peri-menopausal women. Hence I will restrict the time since last menstruation to 
#' mean(length of mentrual cycle) + 2x SD(length of mentrual cycle) = 26.8197 + 2 x 7.34437 = 41.50844
#' Values obtained from [UKB showcase](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=3710) 
#' 
#' Finally, I also want to include only younger women, and restrict the premenopausal group to women younger than 60. Similarily, I want the post-menopausal women to be older, so I can be sure that it is menopause and not some surgical issue. 
#' 

myTab[sex==1, group := "men"]
myTab[menopause==1, group := "post"]
myTab[menopause==0, group := "pre"]
myTab = myTab[!is.na(group)]

myTab[,table(group,gSex)]
myTab = myTab[sex == gSex,]

filt = myTab$group == "pre" & (myTab$lastMenst>41.50844 | myTab$lastMenst<0)
myTab = myTab[!filt,]
filt = myTab$group == "pre" & myTab$age>=60
myTab = myTab[!filt,]

filt = myTab$group == "post" & myTab$age<=50
myTab = myTab[!filt,]

#' ## Regression models ####
#' ***
myTab[smoking == -3, smoking := NA]
myTab[, BMI := log(BMI)]
myTab = myTab[!is.na(BMI) & !is.na(smoking), ]
myTab[is.na(gAncestry) & ancestry != 1001, gAncestry := 0]
myTab[is.na(gAncestry) & ancestry == 1001, gAncestry := 1]

mod1 = lm(PCSK9 ~ smoking + age + lipidMeds + BMI, data = myTab, subset = sex==1)
summary(mod1)
mod2 = lm(PCSK9 ~ smoking + age + lipidMeds + BMI, data = myTab, subset = sex==0)
summary(mod2)

mod3 = lm(PCSK9 ~ smoking + age + lipidMeds + BMI, data = myTab, subset = menopause==0)
summary(mod3)
mod3b = lm(PCSK9 ~ smoking + age + lipidMeds + BMI + lastMenst, data = myTab, subset = menopause==0)
summary(mod3b)

mod4 = lm(PCSK9 ~ smoking + age + lipidMeds + BMI, data = myTab, subset = menopause==1)
summary(mod4)
mod5 = lm(PCSK9 ~ age* menopause + lipidMeds + BMI, data = myTab, subset = sex==0)
summary(mod5)

mod6 = lm(PCSK9 ~ (smoking + age + lipidMeds + BMI)*group + gAncestry, data = myTab)
summary(mod6)

mod7 = lm(PCSK9 ~ (smoking + age + lipidMeds + BMI)*menopause, data = myTab, subset = sex==0)
summary(mod7)

#' **Summary**
#' 
#' - smoking hardly relevant in this data set (weak positive effect throughout)
#' - age has negative effect in men and post-menopausal women, but a positive effect in pre-menopausal women
#' - lipid lowering medication has a positive effect, stronger effect size in men compared to women, and similar in pre- and post-menopausal women
#' - BMI has a positive effect, increasing in size for men - post-menopausal - pre-menopausal women
#' 
#' **Conclusion**
#' 
#' I will keep using the same model as in my previous meta-GWAS, but will stratify for men, post- and pre-menopausal women. 
#' 
#' ## Plots
#' 
#' - boxplot per group
#' - boxplot per ancestry
#' - boxplot per lipid lowering medication and group
#' - scatterplot per BMI and group
#' 
myTab[,strata := paste(lipidMeds,group,sep="_")]
myTab[,lipidMeds2 := "yes"]
myTab[lipidMeds==0,lipidMeds2 := "no"]
myTab[,smoking2 := "2 - current"]
myTab[smoking==1,smoking2 := "1 - previous"]
myTab[smoking==0,smoking2 := "0 - never"]
myTab[,group2 := group]
myTab[group!="men",group2 := paste0(group,"-menopausal women")]
myTab[ancestry< -2,ancestry2 := "7 - prefer not to answer"]
myTab[ancestry< 0 & ancestry >-2,ancestry2 := "8 - do not know"]
myTab[(ancestry< 1005 & ancestry > 10) | ancestry == 1,ancestry2 := "1 - White"]
myTab[(ancestry< 2005 & ancestry > 1005) | ancestry == 2,ancestry2 := "2 - Mixed"]
myTab[(ancestry< 3005 & ancestry > 2005) | ancestry == 3,ancestry2 := "3 - Asian (British)"]
myTab[(ancestry< 4005 & ancestry > 3005) | ancestry == 4,ancestry2 := "4 - Black (British)"]
myTab[ancestry == 5,ancestry2 := "5 - Chinese"]
myTab[ancestry == 6,ancestry2 := "6 - Other"]
myTab[,ancestry3 := ancestry]
myTab[ancestry %in% c(1,2,3,4,5,6),ancestry3 := ancestry * 1000]
myTab[,min(age)]
myTab[,max(age)]
myTab[age<=45, age2 := "40-45"]
myTab[age<=50 & age>45, age2 := "46-50"]
myTab[age<=55 & age>50, age2 := "51-55"]
myTab[age<=60 & age>55, age2 := "56-60"]
myTab[age<=65 & age>60, age2 := "61-65"]
myTab[age<=70 & age>65, age2 := "66-70"]


p1 = ggplot(data = myTab, aes(x = lipidMeds2, y = PCSK9, fill = group)) + 
  facet_wrap(~group2) + 
  geom_boxplot() +
  labs(x="lipid lowering medication",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels by lipid-lowering medication and sex-/age-group")) +
  scale_fill_manual(values = c("steelblue","lightgreen","darkorange"),
                    labels = c("men", "post-menopausal women","premenopausal women"))+
  #theme_classic() +
  theme(legend.position = "none")
p1

p2.1 = ggplot(data = myTab[group == "men"], aes(x = BMI, y = PCSK9, color = group)) + 
  facet_wrap(~lipidMeds2) + 
  geom_point(color = "steelblue") +
  stat_smooth(method='lm', color="darkblue") +
  labs(x="BMI (log-transformed)",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels and BMI by lipid-lowering medication in men")) +
  theme(legend.position = "none")
p2.1

p2.2 = ggplot(data = myTab[group == "post"], aes(x = BMI, y = PCSK9, color = group)) + 
  facet_wrap(~lipidMeds2) + 
  geom_point(color = "lightgreen") +
  stat_smooth(method='lm', color="darkgreen") +
  labs(x="BMI (log-transformed)",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels and BMI by lipid-lowering medication in post-menopausal women")) +
  theme(legend.position = "none")
p2.2

p2.3 = ggplot(data = myTab[group == "pre"], aes(x = BMI, y = PCSK9, color = group)) + 
  facet_wrap(~lipidMeds2) + 
  geom_point(color = "darkorange") +
  stat_smooth(method='lm', color="darkred") +
  labs(x="BMI (log-transformed)",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels and BMI by lipid-lowering medication in pre-menopausal women")) +
  theme(legend.position = "none")
p2.3


p3.1 = ggplot(data = myTab, aes(x = as.factor(ancestry2), y = PCSK9, fill = as.factor(ancestry3))) + 
  geom_boxplot() +
  labs(x="Ethnic background",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels by ethnic background")) +
  #theme_classic() +
  theme(legend.position = "none")
p3.1

p3.2 = ggplot(data = myTab, aes(x = as.factor(ancestry2), y = PCSK9, fill = as.factor(ancestry2))) + 
  geom_boxplot() +
  labs(x="Ethnic background",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels by ethnic background")) +
  #theme_classic() +
  theme(legend.position = "none")
p3.2

p4 = ggplot(data = myTab, aes(x = smoking2, y = PCSK9, fill = group)) + 
  facet_wrap(~group2) + 
  geom_boxplot() +
  labs(x="current smoking",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels by current smoking and sex-/age-group")) +
  scale_fill_manual(values = c("steelblue","lightgreen","darkorange"),
                    labels = c("men", "post-menopausal women","premenopausal women"))+
  #theme_classic() +
  theme(legend.position = "none")
p4

p5.1 = ggplot(data = myTab[group == "men"], aes(x = age2, y = PCSK9, fill = lipidMeds2)) + 
  #facet_wrap(~lipidMeds2) + 
  geom_boxplot() +
  labs(x="age (in years)",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels and age by lipid lowering medication in men")) +
  scale_fill_manual(values = c("steelblue","darkblue"),
                    labels = c("no", "yes"))+
  theme(legend.position = "none")
p5.1

p5.2 = ggplot(data = myTab[group == "post"], aes(x = age2, y = PCSK9, fill = lipidMeds2)) + 
  #facet_wrap(~lipidMeds2) + 
  geom_boxplot() +
  labs(x="age (in years)",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels and age by lipid lowering medication in post-menopausal women")) +
  scale_fill_manual(values = c("lightgreen","darkgreen"),
                    labels = c("no", "yes"))+
  theme(legend.position = "none")
p5.2

p5.3 = ggplot(data = myTab[group == "pre"], aes(x = age2, y = PCSK9, fill = lipidMeds2)) + 
  #facet_wrap(~lipidMeds2) + 
  geom_boxplot() +
  labs(x="age (in years)",
       y="PCSK9 (in NPX)", 
       title = paste0("PCSK9 levels and age by lipid lowering medication in pre-menopausal women")) +
  scale_fill_manual(values = c("darkorange","darkred"),
                    labels = c("no", "yes"))+
  theme(legend.position = "none")
p5.3

ggsave(plot = p1, filename = paste0("../results/_figures/Boxplot_Group_LipLowMed.png"), 
       height = 7, width = 14)
ggsave(plot = p2.1, filename = paste0("../results/_figures/Scatterplot_Men_BMI.png"), 
       height = 7, width = 14)
ggsave(plot = p2.2, filename = paste0("../results/_figures/Scatterplot_Post_BMI.png"), 
       height = 7, width = 14)
ggsave(plot = p2.3, filename = paste0("../results/_figures/Scatterplot_Pre_BMI.png"), 
       height = 7, width = 14)
ggsave(plot = p3.1, filename = paste0("../results/_figures/Boxplot_Ancestry.png"), 
       height = 7, width = 14)
ggsave(plot = p3.2, filename = paste0("../results/_figures/Boxplot_AncestryGrouped.png"), 
       height = 7, width = 14)
ggsave(plot = p4, filename = paste0("../results/_figures/Boxplot_Group_Smoking.png"), 
       height = 7, width = 14)
ggsave(plot = p5.1, filename = paste0("../results/_figures/Boxplot_Men_Age.png"), 
       height = 7, width = 14)
ggsave(plot = p5.2, filename = paste0("../results/_figures/Boxplot_Post_Age.png"), 
       height = 7, width = 14)
ggsave(plot = p5.3, filename = paste0("../results/_figures/Boxplot_Pre_Age.png"), 
       height = 7, width = 14)

#' # Create input for regenie ####
#' ***
#' ## Filter for EUR ancestry
#' 
save(myTab, file = paste0(data_QC,"/01_Prep_01_ukb_PCSK9_filtered.RData"))

myTab = myTab[ancestry2 == "1 - White",]
save(myTab, file = paste0(data_QC,"/01_Prep_01_ukb_PCSK9_filtered_EUR.RData"))

#' ## Sample inclusion 
#' 
#' **Sample inclusion/exclusion file format**: No header. Each line starts with individual FID IID. Space/tab separated.
#' 
#' Needed for the generation of the bed/bim/fam files for regenie step 1
#' 
myTab1 = copy(myTab)
myTab1[,IID := ID]
setnames(myTab1,"ID","FID")
myTab1 = myTab1[,c(1,35)]
write.table(myTab1,
            file = paste0(data_QC,"/01_Prep_01_ukb_SampleList_EUR.txt"), 
            col.names = F, row.names = F, quote = F)

#' ## Covariate file 
#' 
#' **Covariate file format**: Line 1 : Header with FID, IID and C covariate names. Followed by lines of C+2 values. Space/tab separated. Each line contains individual FID and IID followed by C covariate values.
#' 
#' Same file for regenie step 1 and step 2.
#' 
#' This includes: age, smoking, BMI, gArray, gPC1 - 10
#' 
myNames = c("ID", "age", "smoking", "BMI", "gArray", "gPC_1","gPC_2","gPC_3","gPC_4","gPC_5","gPC_6","gPC_7","gPC_8","gPC_9", "gPC_10") 

myTab2 = copy(myTab)
colsOut<-setdiff(colnames(myTab2),myNames)
myTab2[,get("colsOut"):=NULL]
setcolorder(myTab2,myNames)
dim(myTab2)
myTab2[,gArray := gsub("Axiom","1",gArray)]
myTab2[,gArray := gsub("BiLEVE","0",gArray)]
myTab2 = cbind(myTab2$ID,myTab2)
names(myTab2)[1:2] = c("FID","IID")

write.table(myTab2,
            file = paste0(data_QC,"/01_Prep_01_ukb_covariates_EUR.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")

#' ## Phenotype file
#' 
#' **Phenotype file format**: Line 1 : Header with FID, IID and P phenotypes names. Followed by lines of P+2 values. Space/tab separated. Each line contains individual FID and IID followed by P phenotype values (for binary traits, must be coded as 0=control, 1=case, NA=missing unless using --1).
#' 
#' 6 files with one phenotype each for regenie step 1; and 1 file with all six phenotypes for regenie step 2. Reason: step 1 will mean-impute missing phenotypes, but I don't want that for my subgroups as they are independent (no overlaps).  
#' 
myTab3 = copy(myTab)
myTab3[group == "men" & lipidMeds==0,PCSK9_men_free := PCSK9]
myTab3[group == "men" & lipidMeds==1,PCSK9_men_treated := PCSK9]
myTab3[group == "post" & lipidMeds==0,PCSK9_post_free := PCSK9]
myTab3[group == "post" & lipidMeds==1,PCSK9_post_treated := PCSK9]
myTab3[group == "pre" & lipidMeds==0,PCSK9_pre_free := PCSK9]
myTab3[group == "pre" & lipidMeds==1,PCSK9_pre_treated := PCSK9]

myTab3 = cbind(myTab3$ID,myTab3$ID,myTab3[,c(35:40)])
names(myTab3)[1:2] = c("FID","IID")

write.table(myTab3,
            file = paste0(data_QC,"/01_Prep_01_ukb_phenotypes_EUR_step2_PCSK9.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab3[,c(1,2,3)],
            file = paste0(data_QC,"/01_Prep_01_ukb_phenotypes_EUR_step1_PCSK9_1.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab3[,c(1,2,4)],
            file = paste0(data_QC,"/01_Prep_01_ukb_phenotypes_EUR_step1_PCSK9_2.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab3[,c(1,2,5)],
            file = paste0(data_QC,"/01_Prep_01_ukb_phenotypes_EUR_step1_PCSK9_3.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab3[,c(1,2,6)],
            file = paste0(data_QC,"/01_Prep_01_ukb_phenotypes_EUR_step1_PCSK9_4.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab3[,c(1,2,7)],
            file = paste0(data_QC,"/01_Prep_01_ukb_phenotypes_EUR_step1_PCSK9_5.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")
write.table(myTab3[,c(1,2,8)],
            file = paste0(data_QC,"/01_Prep_01_ukb_phenotypes_EUR_step1_PCSK9_6.txt"), 
            col.names = T, row.names = F, quote = F,sep = "\t")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

