metaGWAS_JP = function(data, doGC = T, doREM=T,signDig=6){
  # I expect the following column names in each data set: 
  # beta, SE, pval, sampleSize, EAF, MAF, info, SNP
  # 
  # I expect the data sets are already harmonized (same SNPs, same effect allele)
  # data is a list including the study-specific summary statistics
  
  # Step 1: get meta data
  {
    k = length(data)
    
    ## Step 1a: get unique SNPs
    SNPID = c()
    for (i in 1 : k) {
      SNPID = c(SNPID, as.character(data[[i]][, SNP]))
    }
    
    ## Step 1b: Match data by SNPID
    SNPID = as.factor(unique(SNPID))
    for (i in 1 : k) {
      matched = match(SNPID, data[[i]][, SNP])
      data[[i]] = data[[i]][matched]
    }
    SNPID = as.character(SNPID)
    m = length(SNPID) 
    
    ## Step 1c: get meta data (sample size, EAF, MAF, info)
    n         = c()
    eaf       = c()
    maf       = c()
    infoscore = c()
    for (j in 1 : k) {
      n                        = cbind(n, data[[j]][, sampleSize])
      eaf                      = cbind(eaf, data[[j]][, EAF])
      maf                      = cbind(maf, data[[j]][, MAF])
      infoscore                = cbind(infoscore, data[[j]][, info])
    }
    numberStudies      = apply(!is.na(n), 1, sum)
    totalN             = apply(n, 1, sum, na.rm = TRUE)
    nWeight            = n / totalN
    nWeightedEAF       = apply((nWeight * eaf), 1, sum, na.rm = TRUE)
    nWeightedMAF       = apply((nWeight * maf), 1, sum, na.rm = TRUE)
    nWeightedInfoScore = apply((nWeight * infoscore), 1, sum, na.rm = TRUE)
    
    # Step 1d: merge meta data into a data table
    metaData = data.table(markerID = SNPID,
                          numberStudies = numberStudies,
                          totalN = totalN,
                          nWeightedEAF = nWeightedEAF,
                          nWeightedMAF = nWeightedMAF,
                          nWeightedInfo = nWeightedInfoScore)
    
  }
  
  # Step 2: get effect direction
  direction = rep("", m)
  for (i in 1 : k) {
    dum             = is.na(data[[i]][, beta])
    direction[dum]  = paste0(direction[dum], "?")
    direction[!dum] = paste0(direction[!dum], ifelse((data[[i]][!dum, beta] > 0), "+", "-"))
  }
  metaData$direction = direction
  
  # Step 3: meta-analysis fixed effect
  {
    theta = c()
    w     = c()
    for(j in 1 : k) {
      theta0 = data[[j]][, beta]
      se0    = data[[j]][, SE]
      w0     = 1 / se0^2
      if(doGC==T){
        lambdaGC = median((theta0 / se0)^2, na.rm = TRUE) / 0.456
        if(lambdaGC > 1)   w0 = 1 / (lambdaGC * se0^2)
      }
      theta = cbind(theta, theta0)
      w     = cbind(w, w0)
    }
    thetaFEM = apply((w * theta), 1, sum, na.rm = TRUE) / apply(w, 1, sum, na.rm = TRUE)
    seFEM = sqrt(1 / apply(w, 1, sum, na.rm = TRUE))
    ZFEM = thetaFEM / seFEM
    pFEM = 2 * pnorm(abs(ZFEM),lower.tail = FALSE)
    
    Q       = apply((w * (theta - thetaFEM)^2), 1, sum, na.rm = TRUE)
    I2      = (Q - (2 - 1)) / Q
    dum     = (I2 < 0) | is.na(I2)
    I2[dum] = 0
    
  }
  
  # Step 4: meta-analysis random effect
  if(doREM == T){
    theta = c()
    w     = c()
    for(j in 1 : k) {
      theta0 = data[[j]][, beta]
      se0    = data[[j]][, SE]
      w0     = 1 / se0^2
      theta = cbind(theta, theta0)
      w     = cbind(w, w0)
    }
    thetaFEM2 = apply((w * theta), 1, sum, na.rm = TRUE) / apply(w, 1, sum, na.rm = TRUE)
    dum = numberStudies == 1
    if (sum(dum) > 0) {
      # Function to extract non-NA values of a vector
      extractNonNA  = function(vector) { return (vector[!is.na(vector)]) } 
      # Application of function "extractNonNA" to extract the single theta of each line of the theta matrix 
      # ("numberStudies == 1", i.e. one theta per line, otherwise "NA")
      thetaFEM2[dum] = apply(theta[dum, ], 1, extractNonNA) 
    }
    Q2         = apply((w * (theta - thetaFEM2)^2), 1, sum, na.rm = TRUE)
    I2_REM     = (Q2 - (numberStudies - 1)) / Q2
    dum        = (I2_REM < 0) | is.na(I2_REM)
    I2_REM[dum]= 0
    
    tau2      = (Q2 - (numberStudies - 1)) / (apply(w, 1, sum, na.rm = TRUE) - (apply((w^2), 1, sum, na.rm = TRUE) / apply(w, 1, sum, na.rm = TRUE)))
    dum       = (Q2 < (numberStudies - 1)) | is.na(tau2)
    tau2[dum] = 0
    
    wREM = matrix(nrow = m, ncol = k) 
    for (j in 1 : k) {
      for (i in 1 : m) {
        wREM[i, j] = 1 / (1 / w[i, j] + tau2[i])
      }
    }
    
    thetaREM = apply((wREM * theta), 1, sum, na.rm = TRUE) / apply(wREM, 1, sum, na.rm = TRUE)
    seREM    = sqrt(1 / apply(wREM, 1, sum, na.rm = TRUE))
    ZREM = thetaREM / seREM
    pREM = 2 * pnorm(abs(ZREM), lower.tail = FALSE)
    
  }
  
  # Step 5: merge meta data and meta-analysis results
  metaData$FEM_beta = signif(thetaFEM, digits = signDig)
  metaData$FEM_SE = signif(seFEM, digits = signDig)
  metaData$FEM_pval = signif(pFEM, digits = signDig)
  if(doREM == T){
    metaData$REM_beta = signif(thetaREM, digits = signDig)
    metaData$REM_SE = signif(seREM, digits = signDig)
    metaData$REM_pval = signif(pREM, digits = signDig)
  }
  metaData$CochransQ  = signif(Q, digits = signDig)
  metaData$pCochransQ = signif(pchisq(Q, df = (numberStudies - 1), lower.tail = FALSE), digits = signDig)
  metaData$I2         = signif(I2, digits = signDig)
  
  # Step 6: return result
  return(metaData)
  
}
