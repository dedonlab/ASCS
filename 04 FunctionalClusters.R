############################################################
############################################################
# Functional clusters
# 
# Author: Farzan Ghanegolmohammadi
#
############################################################
############################################################
# Description:
# Hierarchical cluster analysis (HCA) of the obtained functional profile
############################################################
setwd("PATH")
# Prerequisite
if("BiocManager" %in% rownames(installed.packages())){
  library(BiocManager)
} else {install.packages("BiocManager")}

if("qvalue" %in% rownames(installed.packages())){
  library(qvalue)
} else {BiocManager::install("qvalue"); library(qvalue)}
############################################################
## Loading slimmed GO data
load("Data/GO/GOData.rdata")
## Significant level
FDR <- .05
############################################################
# Hierarchical cluster analysis based on functional profile
HCA <- NULL
HCA$hclust <- hclust(dist(GO$ReducedGOmat, method = "binary"), method = "complete")

## Cutting the tree at h = 0.0999
temp <- sort(cutree(HCA$hclust, h = 0.999))
test <- tapply(temp, factor(temp), names)
names(test) <- sprintf("Group%03i", as.numeric(names(test)))
HCA$CutTreeAll <- list(Nodes = temp, Cut = test)
rm(test,temp)

## Funding groups with >2 members
HCA$CutTree_3 <- HCA$CutTreeAll$Cut[sapply(HCA$CutTreeAll$Cut, length) > 2]

save(HCA, file = "Data/GO/HCA.rdata")
######################################
# GO enrichment of the obtained functional groups
#### ORF names from functional grouping given the threshold (>2)
ORF_Names <- unlist(HCA$CutTree_3)

#### GO Boolean matrix with ORF names from functional grouping given the threshold (>2)
GO_Enrich <- GO$ReducedGOmat[ORF_Names,]

Fisher <- NULL
Fisher$PVal <- matrix(NA, nrow = length(HCA$CutTree_3), ncol = ncol(GO_Enrich),
                      dimnames = list(names(HCA$CutTree_3), colnames(GO_Enrich)))
pb <- txtProgressBar(min = 0, max = length(HCA$CutTree_3), style = 3)

for (j in names(HCA$CutTree_3)) {
  GroupedORFs <- vector(mode = "logical", length = length(ORF_Names))
  names(GroupedORFs) <- ORF_Names
  GroupedORFs[HCA$CutTree_3[[j]]] <- TRUE
  if(sum(GroupedORFs) != length(HCA$CutTree_3[[j]])) stop()

  for (i in colnames(GO_Enrich)) {
  temp <- data.frame(y = factor(as.vector(GO_Enrich[,i]), levels = c("TRUE","FALSE")),
                     x = factor(GroupedORFs, levels = c("TRUE","FALSE")))
  # Fisher's exact test
  Fisher$PVal[j,i] <- fisher.test(table(temp), alternative = "greater")$p.value

  rm(temp)
  }
  rm(GroupedORFs)
  setTxtProgressBar(pb, which(names(HCA$CutTree_3) == i))
}
rm(i,j,pb)

Fisher$qVal <- qvalue::qvalue(Fisher$PVal)$qvalues
save(Fisher, file = "Data/GO/FisherTest.rdata")