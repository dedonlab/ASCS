############################################################
############################################################
# Mapping codon-bias ORFs onto the functional network
# 
# Author: Farzan Ghanegolmohammadi
#
############################################################
############################################################
# Description:
## Visualization of codon-bias analogy as a network:
## Pairwise similarities between codon-bias ORFs were calculated as
## Pearson’s correlation coefficients (r) using cCVs.
############################################################
############################################################
setwd("PATH")
# Prerequisite
if("qgraph" %in% rownames(installed.packages())){
  library(qgraph)
} else {install.packages("qgraph"); library(qgraph)}

if("BiocManager" %in% rownames(installed.packages())){
  library(BiocManager)
} else {install.packages("BiocManager")}

if("qvalue" %in% rownames(installed.packages())){
  library(qvalue)
} else {BiocManager::install("qvalue"); library(qvalue)}

## An function for estimating P values for correlation r values
PR <- function(r, n, lower=FALSE) pt((r * sqrt((n - 2) / (1 - r^2))), df=n - 2, lower.tail=lower)
######################################
## Loading CCA results data
load("Data/CCA/CCA.rdata")
## Loading pairwise CCA
load("Data/CCA/CCA_Pair.rdata")
## Loading HCA results (i.e., functional groups)
load("Data/HCA/HCA.rdata")
## Significance threshold
FDR <- .05
FDR_2 <- .001
######################################
# Construction of morphological network by single edge between gene groups
n <- 0.95 * sum(CCA$cor["q",] < FDR)
n <- ifelse(n %% 1 > .5, ceiling(n), floor(n))

## Correlation analysis between gene-pairs based on pairwise CCA results
Pearson_r <- cor(t(CCA_Pair$Data$CVs[unlist(HCA$CutTree_3),1:n,"Codon"]))
## Estimating P value (two-sided, t-test)
Pearson_p <- PR(abs(Pearson_r), ncol(CCA_Pair$Data$CVs)) * 2
## Estimating q values
Pearson_q <- qvalue(Pearson_p)$qvalues
diag(Pearson_q) <- NA
######################################
## within group ("wth"), between group ("btw")
temp <- array(NA, dim=c(length(unlist(HCA$CutTree_3)), length(unlist(HCA$CutTree_3)), 6))
dimnames(temp) <- list(unlist(HCA$CutTree_3), unlist(HCA$CutTree_3), c("wth", "btw", "r","p","qVal","net"))

### Within group
temp[unlist(HCA$CutTree_3),unlist(HCA$CutTree_3),"wth"] <- matrix(0, nrow = length(unlist(HCA$CutTree_3)),
                                                                  ncol = length(unlist(HCA$CutTree_3)),
                                                                  dimnames=list(unlist(HCA$CutTree_3),unlist(HCA$CutTree_3)))

for(i in names(HCA$CutTree_3)){
  temp[HCA$CutTree_3[[i]],HCA$CutTree_3[[i]],"wth"] <- as.numeric(substr(i,6,nchar(i)))
}
rm(i)
### Pearson's r
temp[rownames(Pearson_r),colnames(Pearson_r), "r"] <- Pearson_r
### Pearson's p
temp[rownames(Pearson_q),colnames(Pearson_q), "p"] <- Pearson_p
### Pearson's q
temp[rownames(Pearson_q),colnames(Pearson_q), "qVal"] <- Pearson_q

### Between groups
temp[,,"btw"] <- 0
for(i in 1:nrow(CCA_Pair$Data$Pairs)) {
  if(CCA_Pair$CV1q[CCA_Pair$Data$Pairs[i,"x"], CCA_Pair$Data$Pairs[i,"y"]] < FDR_2) {
    temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]], HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],"btw"] <- i
    temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]], HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],"btw"] <- i

    r <- temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],"r"]
    q <- temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],"qVal"]
    if(sum(q < FDR, na.rm = T) > 0){
      temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],"net"][q < FDR] <- r[q < FDR]
      rm(r,q)
    }else{
      test <- temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],"p"]
      temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],"net"][which(test == min(test), arr.ind=T)] <- temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],"r"][which(test == min(test), arr.ind=T)]
      rm(test)
    }

    r <- temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],"r"]
    q <- temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],"qVal"]
    if(sum(q < FDR, na.rm = T) > 0){
      temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],"net"][q < FDR] <- r[q < FDR]
      rm(r,q)
    }else{
      test <- temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],"p"]
      temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],"net"][which(test == min(test), arr.ind=T)] <- temp[HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"x"]]],HCA$CutTree_3[[CCA_Pair$Data$Pairs[i,"y"]]],"r"][which(test == min(test), arr.ind=T)]
      rm(test)
    }
  }
}
rm(i)

temp[,,"net"][temp[,,"wth"] > 0] <- temp[,,"r"][temp[,,"wth"] > 0]
diag(temp[,,"net"]) <- NA
######################################
Network <- NULL
Network$cormat <- temp
Network$scor <- qgraph(temp[,,"net"], layout = "spring", DoNotPlot = TRUE)
rownames(Network$scor$layout) <- colnames(Network$scor$Arguments$input)

colramp <- colorRamp(c("#0000FF", "#FFFFFF", "#FF0000"))
temp <- data.frame(from = Network$scor$Edgelist$from,
                   to = Network$scor$Edgelist$to,
                   r = Network$scor$Edgelist$weight)

test <- (temp$r + 1)/2
test <- rgb(colramp(test), max = 255)
temp$col <- test
rm(test, colramp)
temp <- temp[order(abs(temp$r)),]

### Adding ORF names to the network
for (i in 1:nrow(Network$scor$layout)) {
  temp$from[which(temp$from == i)] <- rownames(Network$scor$layout)[i]
  temp$to[which(temp$to == i)] <- rownames(Network$scor$layout)[i]
}
rm(i)

## Coordinates of each node
temp$x0 <- Network$scor$layout[temp$from, 1]
temp$y0 <- Network$scor$layout[temp$from, 2]
temp$x1 <- Network$scor$layout[temp$to, 1]
temp$y1 <- Network$scor$layout[temp$to, 2]
Network$scor$Edgelist$color <- temp
######################################
# Visualization
plot(0, 0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), ann = FALSE, axes = FALSE)
segments(x0 = Network$scor$Edgelist$color$x0, y0 = Network$scor$Edgelist$color$y0,
         x1 = Network$scor$Edgelist$color$x1, y1 = Network$scor$Edgelist$color$y1, lwd = .5)
points(Network$scor$layout, pch = 16, col = FinalCol_n, cex = 1, lwd = .1)
points(Network$scor$layout, pch = 1, col = "#00000070", cex = 1, lwd = .1)