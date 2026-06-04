#' ---
#' title: "Get descriptive data"
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
#' I want the classic Table 1 information for all 2 phenotypes x 7 setting: 
#' 
#' - men without statin treatment
#' - men with statin treatment
#' - women with statin treatment
#' - women without statin treatment
#' - post-menopausal women with statin treatment
#' - post-menopausal women without statin treatment
#' - pre-menopausal women without statin treatment
#' 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")

#' # Load UKB Data ####
#' ***
covar_PCSK9 = fread(paste0(data_QC,"/02_regenie_input/01_Prep_ukb_covariates_EUR_PCSK9.txt"))
covar_LDLC = fread(paste0(data_QC,"/02_regenie_input/01_Prep_ukb_covariates_EUR_LDLC.txt"))

phenos_PCSK9 = fread(paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step2_PCSK9.txt"))
phenos_LDLC = fread(paste0(data_QC,"/02_regenie_input/01_Prep_ukb_phenotypes_EUR_step2_LDLC.txt"))

#' Check sample overlap
#' 
table(phenos_LDLC$FID %in% phenos_PCSK9$IID)
table(phenos_PCSK9$FID %in% phenos_LDLC$IID)

#' All good!
#' 
#' # Check LDLC data ####
#' ***
#' 
myPhenos = names(phenos_LDLC)[3:9]
dumTab = foreach(i = 1:length(myPhenos))%do%{
  #i=1
  message("Working on phenotype ",myPhenos[i])
  
  # Split trait
  phenotype = gsub("_.*","",myPhenos[i])
  statin = gsub(".*_","",myPhenos[i])
  sex = gsub(paste0(phenotype,"_"),"",myPhenos[i])
  sex = gsub(paste0("_",statin),"",sex)
  
  # Get phenotype information
  phenos_LDLC[,LDLC := get(myPhenos[i])]
  mean_LDLC = phenos_LDLC[,mean(LDLC,na.rm=T)]
  sd_LDLC = phenos_LDLC[,sd(LDLC,na.rm=T)]
  n_LDLC = dim(phenos_LDLC[!is.na(LDLC)])[1]
  
  res = data.table(phenotype = phenotype, 
                   sex = sex, 
                   statin = statin, 
                   trait = myPhenos[i],
                   n_trait = n_LDLC,
                   mean_trait = mean_LDLC,
                   sd_trait = sd_LDLC)
  res
  
  # get covariables information
  covar = copy(covar_LDLC)
  covar = covar[IID %in% phenos_LDLC[!is.na(LDLC),IID]]
  
  res[,mean_age := covar[,mean(age)]]
  res[,sd_age := covar[,sd(age)]]
  res[,mean_logBMI := covar[,mean(BMI)]]
  res[,sd_logBMI := covar[,sd(BMI)]]
  res[,n_smoke_abs := dim(covar[smoking == 1,])[1]]
  res[,n_smoke_rel := 100 * n_smoke_abs / n_LDLC]
  res 
}
tab_LDLC = rbindlist(dumTab)
kable(tab_LDLC,digits = 2)

#' # Check PCSK9 data ####
#' ***
#' 
myPhenos = names(phenos_PCSK9)[3:9]
dumTab = foreach(i = 1:length(myPhenos))%do%{
  #i=1
  message("Working on phenotype ",myPhenos[i])
  
  # Split trait
  phenotype = gsub("_.*","",myPhenos[i])
  statin = gsub(".*_","",myPhenos[i])
  sex = gsub(paste0(phenotype,"_"),"",myPhenos[i])
  sex = gsub(paste0("_",statin),"",sex)
  
  # Get phenotype information
  phenos_PCSK9[,PCSK9 := get(myPhenos[i])]
  mean_PCSK9 = phenos_PCSK9[,mean(PCSK9,na.rm=T)]
  sd_PCSK9 = phenos_PCSK9[,sd(PCSK9,na.rm=T)]
  n_PCSK9 = dim(phenos_PCSK9[!is.na(PCSK9)])[1]
  
  res = data.table(phenotype = phenotype, 
                   sex = sex, 
                   statin = statin, 
                   trait = myPhenos[i],
                   n_trait = n_PCSK9,
                   mean_trait = mean_PCSK9,
                   sd_trait = sd_PCSK9)
  res
  
  # get covariables information
  covar = copy(covar_PCSK9)
  covar = covar[IID %in% phenos_PCSK9[!is.na(PCSK9),IID]]
  
  res[,mean_age := covar[,mean(age)]]
  res[,sd_age := covar[,sd(age)]]
  res[,mean_logBMI := covar[,mean(BMI)]]
  res[,sd_logBMI := covar[,sd(BMI)]]
  res[,n_smoke_abs := dim(covar[smoking == 1,])[1]]
  res[,n_smoke_rel := 100 * n_smoke_abs / n_LDLC]
  res[,mean_timeDif := covar[,mean(PPP_timeDif)]]
  res[,sd_timeDif := covar[,sd(PPP_timeDif)]]
  
  res 
}
tab_PCSK9 = rbindlist(dumTab)
kable(tab_PCSK9,digits = 2)

#' # Save tables ####
#' ***
tab = rbind(tab_PCSK9,tab_LDLC,fill=T)
tosave4 = data.table(data = c("tab", "tab_PCSK9", "tab_LDLC"), 
                     SheetNames = c("Merged","PCSK9","LDLC"))
excel_fn = paste0("../results/05_descriptiveTables.xlsx")
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(tab, tab_PCSK9, tab_LDLC, 
     file = paste0("../results/05_descriptiveTables.RData"))

write.table(tab, file = paste0("../results/05_descriptiveTables.txt"),
            col.names = T,row.names = F,quote = F,sep="\t")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

