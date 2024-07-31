############################################################
############################################################
# Functional profile
# 
# Author: Farzan Ghanegolmohammadi
#
############################################################
############################################################
# Description:
## Preparing functional profile:
## The basic version of the GO (go-basic.obo) was downloded form Gene Ontology Consortium (http://geneontology.org/docs/download-ontology/).
## Gene annotations were downloaded from Saccharomyces Genome Database (SGD; https://www.yeastgenome.org/).
## A Boolean matrix of GO terms was generated in which if an ORF was annotated by a GO, the value was TRUE; otherwise FALSE.
## Next, a GO slimmer process was done in five steps:
### 1) Removing global GO terms (i.e., >200 repetitions)
### 2) Removing GO terms with identical sets of annotated ORFs
### 3) Removing non-codon-bias ORFs
### 4) Removing unique GO terms (i.e., <2 repetitions)
### 5) Removing ORFs with no annotations
############################################################
############################################################
setwd("PATH")
# Prerequisite
source("Functions/GO_Slimer.r")
############################################################
## Codon information
CodonInfo <- read.csv("Data/CodonInfo.csv", header = TRUE, row.names = 1)
dim(CodonInfo)
## 59 codons
C59 <- CodonInfo$Codon[!CodonInfo$Codon %in% c("ATG", "TGG", "TAA","TAG","TGA")]
## BI q values
load("Data/GLLM/BIqValues.rdata")
## GO data
GO <- NULL
### GO Boolean matrix
GO$GOmat <- read.table("Data/GO/GOmatrix.tsv", header = TRUE,
                       row.names = 1, sep = "\t", stringsAsFactors = F)
colnames(GO$GOmat) <- chartr(".",":", colnames(GO$GOmat))
GO$GOmat <- as.matrix(GO$GOmat)
### GO terms
GO$terms <- read.csv("Data/GO/GOterms.csv", header = TRUE, row.names = 1)
sapply(GO, dim)

## FDR threshold
FDR <- 0.05
############################################################
## No. of significant ORFs (q < FDR) in each codon
SigCodons <- colnames(BIqValues)[colSums(BIqValues < FDR, na.rm = T) > 0]
all(SigCodons == C59) # All of C59 have at least one significant ORF

## No. of Significant codons (q < FDR) in each ORF
SigORFs <- rownames(BIqValues)[rowSums(BIqValues < FDR, na.rm = T) > 0]
length(SigORFs) # 1776
############################################################
# GO slimmer steps:
## 1) Removing global GO terms (i.e., >200 repetitions)
GO$ReducedGOmat <- GO$GOmat[,colnames(GO$GOmat)[apply(GO$GOmat, MARGIN = 2, sum) < 200]]
## 2) Removing GO terms with identical sets of annotated ORFs
GO_Slim <- goSlimer(GO$ReducedGOmat)
GO$ReducedGOmat <- GO$ReducedGOmat[,names(GO_Slim)]
## 3) Removing non-codon-bias ORFs
GO$ReducedGOmat <- GO$ReducedGOmat[rownames(GO$ReducedGOmat) %in% SigORFs,]
## 4) Removing unique GO terms (i.e., <2 repetitions)
GO$ReducedGOmat <- GO$ReducedGOmat[,colnames(GO$ReducedGOmat)[apply(GO$ReducedGOmat, MARGIN = 2, sum) > 2]]
## 5) Removing ORFs with no annotations
GO$ReducedGOmat <- GO$ReducedGOmat[rowSums(GO$ReducedGOmat) > 0,]

save(GO, file = "Data/GO/GOData.rdata")