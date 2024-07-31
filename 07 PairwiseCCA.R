############################################################
############################################################
# Associations among functional clusters
# 
# Author: Farzan Ghanegolmohammadi
# 
############################################################
############################################################
# Description:
## Pairwise Canonical Correlation Analysis (CCA) was implemented in following steps:
## 1) To avoid over-fitting, Principal Component Analysis (PCA) was applied to
##    codon-bias canonical variables (cCV) scores of members of each cluster.
## 2) CCA was applied to the obtained PC scores of every pair, where the cumulative
##    number of PC scores of every pair is less than the number of significant CVs.
## 3) The obtained heterogeneous canonical variables (hCVs) served as independent
##    components correlating functional clusters.
## The significance of the canonical correlation coefficient of the first hCVs was tested
## using Bartlett’s chi-squared test.
############################################################
############################################################
setwd("PATH")
# Prerequisite
source("Functions/pca.r")
source("Functions/CanCor.r")
############################################################
# ## Loading CCA results data
load("Data/CCA/CCA.rdata")
# ## Loading HCA results (i.e., functional Clusters)
load("Data/GO/HCA.rdata")
## Significant levels
FDR <- .05
FDR_2 <- .001
############################################################
# Pairwise CCA (CCA among functional clusters)
CCA_Pair <- NULL
Codon <- CCA$Codon$score[unlist(HCA$CutTree_3), which(CCA$cor["q",] < FDR)]
GO_res <- CCA$GO$score[unlist(HCA$CutTree_3), which(CCA$cor["q",] < FDR)]
CCA_Pair$Data$CVs <- array(NA, dim = c(nrow(Codon), ncol(Codon), 2))
dimnames(CCA_Pair$Data$CVs) <- list(rownames(Codon), colnames(CCA$cor)[CCA$cor["q",] < FDR], c("Codon","GO"))

CCA_Pair$Data$CVs[rownames(Codon),,"Codon"] <- Codon
CCA_Pair$Data$CVs[rownames(GO_res),,"GO"] <- GO_res
sum(is.na(CCA_Pair$Data$CVs)) == 0 # TRUE
rm(Codon, GO_res)

# No. of employed dimensions
## To avoid over fitting, we use less of no. of the significant canonical variables (CVs) at FDR == FDR
n <- 0.95 * sum(CCA$cor["q",] < FDR)
n <- ifelse(n %% 1 > .5, ceiling(n), floor(n))

### Choosing pairs of Clusters
temp <- NULL
for(i in 1:(length(HCA$CutTree_3)-1)) {
  for(j in (i+1):(length(HCA$CutTree_3))) {
    x <- names(HCA$CutTree_3)[i]
    y <- names(HCA$CutTree_3)[j]
    temp[[paste0(x, sub("Cluster", "vs", y))]] <- c(x, y)
  }
  rm(x,y)
}
rm(i,j)

temp <- data.frame(t(data.frame(temp)))
colnames(temp) <- c("x", "y")
CCA_Pair$Data$Pairs <- temp
rm(temp)
############################################################
# Pairwise CCA
pb <- txtProgressBar(min = 0, max = nrow(CCA_Pair$Data$Pairs), style = 3)

res <- NULL
for(i in rownames(CCA_Pair$Data$Pairs)) {
  temp <- NULL
  temp$i <- CCA_Pair$Data$Pairs[i,"x"]
  temp$j <- CCA_Pair$Data$Pairs[i,"y"]
  temp$x <- t(CCA_Pair$Data$CVs[HCA$CutTree_3[[temp$i]],,"Codon"])
  res[[i]]$x <- temp$x
  temp$y <- t(CCA_Pair$Data$CVs[HCA$CutTree_3[[temp$j]],,"Codon"])
  res[[i]]$y <- temp$y
  if(ncol(temp$x) + ncol(temp$y) > n) {
    if(ncol(temp$x) > n/2 & ncol(temp$y) > n/2) {
      temp$PCA$x <- pca(temp$x)
      temp$PCA$y <- pca(temp$y)
      temp$PCA$CCR$x <- temp$PCA$x$importance["Cumulative Proportion",]
      names(temp$PCA$CCR$x) <- paste0("x", names(temp$PCA$CCR$x))
      temp$PCA$CCR$y <- temp$PCA$y$importance["Cumulative Proportion",]
      names(temp$PCA$CCR$y) <- paste0("y", names(temp$PCA$CCR$y))
      temp$PCA$CCR <- sort(c(temp$PCA$CCR$x, temp$PCA$CCR$y))
      temp$PCA$CCR <- temp$PCA$CCR[1:n]
      temp$PCA$CCR <- list(x = sub("x", "", names(temp$PCA$CCR)[grep("^x", names(temp$PCA$CCR))]),
                           y = sub("y", "", names(temp$PCA$CCR)[grep("^y", names(temp$PCA$CCR))]))
      res[[i]]$PCA <- temp$PCA
      temp$res <- CanCor(temp$PCA$x$x[,temp$PCA$CCR$x, drop = F], temp$PCA$y$x[,temp$PCA$CCR$y, drop = F])
    } else if (ncol(temp$x) < ncol(temp$y)) {
      temp$PCA <- pca(temp$y)
      res[[i]]$PCA$x <- NA
      res[[i]]$PCA$y <- temp$PCA
      res[[i]]$PCA$CCR$x <- NA
      res[[i]]$PCA$CCR$y <- colnames(temp$PCA$x)[1:(n-ncol(temp$x))]
      temp$res <- CanCor(temp$x, temp$PCA$x[,1:(n-ncol(temp$x)), drop = F])
    } else {
      temp$PCA <- pca(temp$x)
      res[[i]]$PCA$x <- temp$PCA
      res[[i]]$PCA$y <- NA
      res[[i]]$PCA$CCR$x <- colnames(temp$PCA$x)[1:(n-ncol(temp$y))]
      res[[i]]$PCA$CCR$y <- NA
      temp$res <- CanCor(temp$PCA$x[,1:(n-ncol(temp$y)), drop = F], temp$y)
    }
  } else {
    res[[i]]$PCA$x <- NA
    res[[i]]$PCA$y <- NA
    res[[i]]$PCA$CCR$x <- NA
    res[[i]]$PCA$CCR$y <- NA
    temp$res <- CanCor(temp$x, temp$y)
  }
  res[[i]]$cca <- temp$res
  setTxtProgressBar(pb, which(rownames(CCA_Pair$Data$Pairs) == i))
}
CCA_Pair$res <- res
rm(i,pb,res)
############################################################
## P value
CCA_Pair$CV1p <- matrix(NA, nrow=length(HCA$CutTree_3), ncol=length(HCA$CutTree_3),
                          dimnames = list(names(HCA$CutTree_3), names(HCA$CutTree_3)))
## q value
CCA_Pair$CV1q <- matrix(NA, nrow=length(HCA$CutTree_3), ncol=length(HCA$CutTree_3),
                        dimnames = list(names(HCA$CutTree_3), names(HCA$CutTree_3)))

for(i in names(CCA_Pair$res)) {
  ## P value
  CCA_Pair$CV1p[CCA_Pair$Data$Pairs[i,"x"], CCA_Pair$Data$Pairs[i,"y"]] <- CCA_Pair$res[[i]]$cca$cor["p","CV1"]
  CCA_Pair$CV1p[CCA_Pair$Data$Pairs[i,"y"], CCA_Pair$Data$Pairs[i,"x"]] <- CCA_Pair$res[[i]]$cca$cor["p","CV1"]

  ## q value
  CCA_Pair$CV1q[CCA_Pair$Data$Pairs[i,"x"], CCA_Pair$Data$Pairs[i,"y"]] <- CCA_Pair$res[[i]]$cca$cor["q","CV1"]
  CCA_Pair$CV1q[CCA_Pair$Data$Pairs[i,"y"], CCA_Pair$Data$Pairs[i,"x"]] <- CCA_Pair$res[[i]]$cca$cor["q","CV1"]
}
rm(i)
############################################################
## P value correction
colSums(CCA_Pair$CV1q < FDR_2, na.rm = T)
barplot(colSums(CCA_Pair$CV1q < FDR_2, na.rm = T), las = 2)
range(colSums(CCA_Pair$CV1q < FDR_2, na.rm=T))

## Saving the results
save(CCA_Pair, file = "Data/CCA/CCA_Pair.rdata")
