############################################################
############################################################
# Canonical Correlation Analysis (CCA)
# 
# Author: Farzan Ghanegolmohammadi
#
############################################################
############################################################
# Description:
## A two-step strategy were used for CCA:
## 1) Principal Component Analysis (PCA):
## 1-1) Codon-bias profile was subjected to PCA.PCs accounting for 95% of variation of
##      the codon-bias profile are referred to as codon-bias principal components (cPCs).
## 1-2) Functional profile was subjected to PCA. PCs accounting for 99% of variation of
##      the functional profile are referred to as GO principal components (gPCs).
## 2) Canonical Correlation Analysis (CCA):
## 2-1) CCA was applied to cPCs and gPCs.
## 2-2) The significance of the canonical variables was tested using Bartlett’s Chi-square test (one-sided, FDR = 0.05)
##      for both codon-bias canonical variables (cCV) and GO canonical variables (gCV).
############################################################
############################################################
setwd("PATH")
# Prerequisite
source("Functions/CCA.r")
############################################################
# Reading data into R
## Codon information
CodonInfo <- read.csv("Data/CodonInfo.csv", header = TRUE, row.names = 1)
## 59 codons
C59 <- CodonInfo$Codon[!CodonInfo$Codon %in% c("ATG", "TGG", "TAA","TAG","TGA")]
## BI Z values
load("Data/GLM/BIZValues.rdata")
## BI q values
load("Data/GLM/BIqValues.rdata")
## Loading GO data
load("Data/GO/GOData.rdata")
dim(GO$GOmat)
## FDR threshold
FDR <- 0.05
############################################################
## Finding significant codons at FDR
SigCodon <- colnames(BIqValues)[colnames(BIqValues < FDR) > 0]
############################################################
# PCA
PCA <- NULL
## 1) PCA on GO data
PCA$GO$PCA <- pca(GO$ReducedGOmat)
PCA$GO$CCR99 <- colnames(PCA$GO$PCA$importance)[PCA$GO$PCA$importance["Cumulative Proportion",] < .99]
# Preparing Z values
PCA$Codon$Data <- BIZValues[rownames(GO$ReducedGOmat),C59]
## PCA on Z values
PCA$Codon$PCA <- pca(PCA$Codon$Data)
PCA$Codon$CCR95 <- colnames(PCA$Codon$PCA $importance)[PCA$Codon$PCA$importance["Cumulative Proportion",] < .95]
############################################################
# CCA
CCA <- CanCor(x = PCA$GO$PCA$x[, PCA$GO$CCR99],
              y = PCA$Codon$PCA$x[, PCA$Codon$CCR95],
              xpfix = "GO", ypfix = "Codon")
# Saving the results
save(CCA, file = "Data/CCA/CCA.rdata")
