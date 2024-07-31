############################################################
############################################################
# Statistical modeling
# 
# Author: Farzan Ghanegolmohammadi
#
############################################################
############################################################
# Description:
## Implementation of Generalized Liner Model (GLM):
## Codon counts were modeled given their synonymous codon(s) using Binomial distribution
############################################################
############################################################
setwd("PATH/")
# Prerequisite

if("gamlss" %in% rownames(installed.packages())){
  library(gamlss)
} else {install.packages("gamlss"); library(gamlss)}

if("BiocManager" %in% rownames(installed.packages())){
  library(BiocManager)
} else {install.packages("BiocManager")}

if("qvalue" %in% rownames(installed.packages())){
  library(qvalue)
} else {BiocManager::install("qvalue"); library(qvalue)}

## Logistic function
logistic <- function(x) 1/(1+exp(-x))

source("Functions/Likelihood.r")
source("Functions/Summary.r")
source("Functions/vcoc-gamlss.r")
############################################################
## Codon information
CodonInfo <- read.csv("Data/CodonInfo.csv", header = TRUE, row.names = 1)
## 59 codons
C59 <- CodonInfo$Codon[!CodonInfo$Codon %in% c("ATG", "TGG", "TAA","TAG","TGA")]
# Reading count data of each codon of yeast (S. cerevisiae strain S288C)
CodonCount_n <- read.csv("Data/CodonCount_n.csv", header = TRUE, row.names = 1)
# Reading count data of codon of yeast (S. cerevisiae strain S288C)
CodonCount_N <- read.csv("Data/CodonCount.csv", header = TRUE, row.names = 1)
## FDR threshold
FDR <- 0.05
############################################################
# Fitting a binomial model
gcont <- gamlss.control(n.cyc = 200L, trace = FALSE)
BIfit <- NULL
for (i in C59) {
  temp <- NULL
  temp <- data.frame(n = CodonCount_n[,i], N = CodonCount_N[,i])
  BIfit[[i]] <- gamlss(cbind(n, N-n) ~ 1, data = na.omit(temp), family = BI, control = gcont)
}
rm(i,temp)
############################################################
# P value estimation
BIPValues <-  matrix(data = NA ,nrow = nrow(CodonCount_n), ncol = ncol(CodonCount_n),
                     dimnames = list(rownames(CodonCount_n), colnames(CodonCount_n)))
pb <- txtProgressBar(min = 0, max = length(C59), style = 3)

for (i in C59) {
  temp <- res <- res_corr <- NULL
  ## Data
  temp <- data.frame(n = CodonCount_n[,i], N = CodonCount_N[,i], row.names = rownames(CodonCount_n))
  res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(temp), c("LowerTail","UpperTail")))
  
  for (k in rownames(temp)) {
    ### When a given AA does not exist in the protein sequence (n == 0 and N == 0)
    if(temp[k,"n"] == 0 & temp[k,"N"] == 0){
      res[k,"LowerTail"] <- NA
      res[k,"UpperTail"] <- NA
      next()
    }
    ### Lower tail
    res[k,"LowerTail"] <- pBI(q = temp[k,"n"], bd = temp[k,"N"], mu = logistic(BIfit[[i]]$mu.coefficients), lower.tail = T)
    ### Upper tail
    ## Correcting for n == 0 or  n == N
    if(temp[k,"n"] == 0 | temp[k,"n"] == temp[k,"N"]){ 
      ### Lower.tail
      res[k,"UpperTail"] <- pBI(q = temp[k,"n"], bd = temp[k,"N"], mu = logistic(BIfit[[i]]$mu.coefficients), lower.tail = F) +
        dBI(x = temp[k,"n"], bd = temp[k,"N"], mu = logistic(BIfit[[i]]$mu.coefficients), log = F)
    }else{ ## when n!= 0 or n != N
      res[k,"UpperTail"] <- pBI(q = temp[k,"n"], bd = temp[k,"N"], mu = logistic(BIfit[[i]]$mu.coefficients), lower.tail = F) 
    }
  }
  ### Two-sided test
  P2Sided <- 2 * (apply(res, MARGIN = 1, min))
  BIPValues[names(P2Sided),i] <- P2Sided
  rm(temp, P2Sided, k, res)
  
  setTxtProgressBar(pb, which(C59 == i))
}
rm(i,pb)
save(BIPValues, file = "Data/GLM/BIPValues.rdata")
###########################################################
# q value estimation
BIPValues[BIPValues > 1] <- 1
BIqValues <- qvalue(BIPValues[,C59])$qvalues
BIqValues[is.na(BIqValues)] <- 1
save(BIqValues, file = "Data/GLM/BIqValues.rdata")
############################################################
## Codon-bias ORFs (q < FDR)
SigORFs <- rownames(BIqValues)[rowSums(BIqValues < FDR, na.rm = T) > 0]
############################################################
# Z value estimation
BIZValues <-  matrix(data = NA ,nrow = nrow(CodonCount_n), ncol = ncol(CodonCount_n),
                     dimnames = list(rownames(CodonCount_n), colnames(CodonCount_n)))
pb <- txtProgressBar(min = 0, max = length(C59), style = 3)

for (i in C59) {
  setTxtProgressBar(pb, which(C59 == i))
  temp <- fitGLM <- NULL
  # Weight for n == 0 or n == N
  w <- ifelse(CodonCount_n[,i] == 0 | CodonCount_n[,i] == CodonCount_N[,i], FALSE, TRUE)
  
  temp <- data.frame(n = CodonCount_n[,i], N = CodonCount_N[,i],
                     x = factor(rownames(CodonCount_n), levels = rownames(CodonCount_n)),
                     w = w, o = rep(BIfit[[i]]$mu.coefficients, nrow(CodonCount_n)))
  # GLM
  fitGLM <- gamlss(cbind(n, N-n) ~ x + offset(o) - 1, data = temp, family = BI, weights = temp$w, control = gcont)
  
  TableANOVA <- summary(fitGLM)
  rownames(TableANOVA) <- sub("x","", rownames(TableANOVA))
  # Saving Z-values
  BIZValues[rownames(TableANOVA),i] <- TableANOVA[,grep("(t|z) value", colnames(TableANOVA))]
  
  # Correcting for n == 0 or n == N
  if(sum(is.na(BIZValues[,i])) > 0){
    for(k in rownames(BIZValues)[is.na(BIZValues[,i])]) {
      n <- qBI(0.5, bd = temp$N[temp$x == k], mu = logistic(BIfit[[i]]$mu.coefficients))
      
      if(temp$n[temp$x == k] == n) {
        BIZValues[k,i] <- 0
      } else if(temp$n[temp$x == k] == 0) {
        if(min(TableANOVA[,grep("(t|z) value", colnames(TableANOVA))]) > 0) {
          BIZValues[k,i] <- 0
        } else {
          BIZValues[k,i] <- min(TableANOVA[,grep("(t|z) value", colnames(TableANOVA))])
        }
      } else {
        if(max(TableANOVA[,grep("(t|z) value", colnames(TableANOVA))]) < 0) {
          BIZValues[k,i] <- 0
        } else {
          BIZValues[k,i] <- max(TableANOVA[,grep("(t|z) value", colnames(TableANOVA))])
        }
      }
    }
  }
}
rm(i,temp,fitGLM,f)
save(BIZValues, file = "Data/GLM/BIZValues.rdata")
