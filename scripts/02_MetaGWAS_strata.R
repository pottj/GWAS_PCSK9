#' ---
#' title: "Run meta-analysis"
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
source("../helperfunction/metaGWAS_JP.R")

#' # Load Summary Statistics ####
#' ***
myFiles = list.files(path = paste0(data_QC,"/regenie/"),pattern = "ukb_step2_PCSK9_chr")
myFiles = myFiles[!grepl(".log",myFiles)]

myTraits = gsub("ukb_step2_PCSK9_chr","",myFiles)
myTraits = gsub(".regenie","",myTraits)
myTraits = gsub("[1234567890]","",myTraits)
myTraits = gsub("_PCSK","PCSK9",myTraits)
myTraits = unique(myTraits)
myTraits

dumTab1 = foreach(i = 1:length(myTraits))%do%{
  #i=1
  myFiles2 = myFiles[grep(myTraits[i],myFiles)]
  
  dumTab2 = foreach(j = 1:length(myFiles2))%do%{
    erg1 = fread(paste0(data_QC,"/regenie/",myFiles2[j]))
    erg1
  }
  erg2 = rbindlist(dumTab2)
  erg2[,phenotype := myTraits[i]]
  setorder(erg2,CHROM,GENPOS)
  
  # Filter rare and low info variants 
  erg2 = erg2[A1FREQ >= 0.01,]
  erg2 = erg2[A1FREQ <= 0.99,]
  erg2 = erg2[INFO >= 0.5,]
  
  # Update column names
  erg2[,TEST := NULL]
  erg2[,CHISQ := NULL]
  erg2[,EXTRA := NULL]
  erg2[,markername := paste(ID, GENPOS, ALLELE0, ALLELE1, sep=":")]
  erg2[!grepl("rs",ID),markername := paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep=":")]
  names(erg2) = c("chr","bp_hg19","rsID", "OA","EA", "EAF","info","nSamples", "beta","SE","logP","phenotype","markername")
  erg2 = erg2[,c(13,3,1,2,4:12)]
  
  # save data in a zipped file
  filenm = paste0(data_QC,"/SummaryStatistics/",myTraits[i],".gz")
  fwrite(erg2, file = filenm)
  
  erg2

}

#' # Meta-analysis ####
#' ***
#' I want to test the following combinations: 
#' 
#' - men: men-free & men-treated
#' - post: post-free & post-treated
#' - pre: pre-free & pre-treated
#' - women: post-free & post-treated & pre-free & pre-treated
#' - women-free: post-free & pre-free
#' - women-treated: post-treated & pre-treated
#' - statin-free: post-free & pre-free & men-free
#' - statin-treated: post-treated & pre-treated & men-treated
#' - overall: men-free & men-treated & post-free & post-treated & pre-free & pre-treated
#' 
ToDoList = data.table(phenotype = c("men","post","pre","women","women_free","women_treated","free","treated","overall"), 
                      dataset1 = myTraits[c( 1, 3, 5,3, 3, 4, 3, 4,1)], 
                      dataset1_NR = c(1,3,5,3,3,4,3,4,1),
                      dataset2 = myTraits[c( 2, 4, 6,4, 5, 6, 5, 6,2)], 
                      dataset2_NR = c(2, 4, 6,4, 5, 6, 5, 6,2),
                      dataset3 = myTraits[c(NA,NA,NA,5,NA,NA, 1, 2,3)],
                      dataset3_NR = c(NA,NA,NA,5,NA,NA, 1, 2,3),
                      dataset4 = myTraits[c(NA,NA,NA,6,NA,NA,NA,NA,4)],
                      dataset4_NR = c(NA,NA,NA,6,NA,NA,NA,NA,4))
ToDoList[,phenotype := paste0("PCSK9_",phenotype)]

for(i in 1:dim(ToDoList)[1]){
  #i=4
  myRow = ToDoList[i,]
  
  # get data sets
  erg1 = dumTab1[[myRow$dataset1_NR]]
  erg2 = dumTab1[[myRow$dataset2_NR]]
  if(!is.na(myRow$dataset3_NR)) erg3 = dumTab1[[myRow$dataset3_NR]]
  if(!is.na(myRow$dataset4_NR)) erg4 = dumTab1[[myRow$dataset4_NR]]
  if(myRow$phenotype=="overall"){
    erg5 = dumTab1[[5]]
    erg6 = dumTab6[[6]]
  }
  # change column names
  setnames(erg1,"nSamples","sampleSize")
  setnames(erg1,"markername","SNP")
  setnames(erg2,"nSamples","sampleSize")
  setnames(erg2,"markername","SNP")
  if(!is.na(myRow$dataset3_NR)){
    setnames(erg3,"nSamples","sampleSize")
    setnames(erg3,"markername","SNP")
  }  
  if(!is.na(myRow$dataset4_NR)){
    setnames(erg4,"nSamples","sampleSize")
    setnames(erg4,"markername","SNP")
  }  
  if(myRow$phenotype=="overall"){
    names(erg5) = names(erg1)
    names(erg6) = names(erg1)
  }
  
  # add MAF
  erg1[,MAF := EAF]
  erg1[EAF>0.5,MAF := 1-EAF]
  erg2[,MAF := EAF]
  erg2[EAF>0.5,MAF := 1-EAF]
  if(!is.na(myRow$dataset3_NR)){
    erg3[,MAF := EAF]
    erg3[EAF>0.5,MAF := 1-EAF]
  }  
  if(!is.na(myRow$dataset4_NR)){
    erg4[,MAF := EAF]
    erg4[EAF>0.5,MAF := 1-EAF]
  }  
  if(myRow$phenotype=="overall"){
    erg5[,MAF := EAF]
    erg5[EAF>0.5,MAF := 1-EAF]
    erg6[,MAF := EAF]
    erg6[EAF>0.5,MAF := 1-EAF]
  }
  
  # merge in one list
  myData1 = list(erg1,erg2)
  if(!is.na(myRow$dataset3_NR) & is.na(myRow$dataset4_NR)) myData1 = list(erg1,erg2,erg3)
  if(!is.na(myRow$dataset4_NR) & myRow$phenotype!= "overall") myData1 = list(erg1,erg2,erg3,erg4)
  if(myRow$phenotype=="overall") myData1 = list(erg1,erg2,erg3,erg4,erg5,erg6)
    
  # run meta-analysis
  erg_meta = metaGWAS_JP(data = myData1,doGC = F,doREM = F,signDig = 6)
  
  # add SNP information
  mySNPs = rbind(erg1[,1:6],erg2[,1:6])
  if(!is.na(myRow$dataset3_NR) & is.na(myRow$dataset4_NR)) mySNPs = rbind(erg1[,1:6],erg2[,1:6],erg3[,1:6])
  if(!is.na(myRow$dataset4_NR) & myRow$phenotype!= "overall") mySNPs = rbind(erg1[,1:6],erg2[,1:6],erg3[,1:6],erg4[,1:6])
  if(myRow$phenotype=="overall") mySNPs = rbind(erg1[,1:6],erg2[,1:6],erg3[,1:6],erg4[,1:6],erg5[,1:6],erg6[,1:6])
  
  mySNPs = mySNPs[!duplicated(SNP),]
  setorder(mySNPs,chr,bp_hg19)
  stopifnot(erg_meta$markerID %in% mySNPs$SNP)
  matched = match(mySNPs$SNP,erg_meta$markerID)
  erg_meta = erg_meta[matched,]
  stopifnot(mySNPs$SNP == erg_meta$markerID)
   
  erg_meta = cbind(mySNPs, erg_meta[,c(4,6,2,3,7:10,13)]) 
  erg_meta[,phenotype := myRow$phenotype]
  setnames(erg_meta,"SNP","markername")
  setnames(erg_meta,"nWeightedEAF","EAF")
  setnames(erg_meta,"nWeightedInfo","info")
  setnames(erg_meta,"numberStudies","nStudies")
  setnames(erg_meta,"totalN","nSamples")
  names(erg_meta) = gsub("FEM_","",names(erg_meta))
  names(erg_meta)
  
  #' save 
  filenm = paste0(data_QC,"/SummaryStatistics/",myRow$phenotype,".gz")
  fwrite(erg_meta, file = filenm)
  
}

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
