############################################################
############################################################
# Canonical Correlation Analysis
# 
# Author: Shinsuke Ohnuki
# 
############################################################
############################################################
## A package for downloading package from bioconductor
if("BiocManager" %in% rownames(installed.packages())){
  library(BiocManager)
} else {install.packages("BiocManager"); library(BiocManager)}
## A package to estimate q value
if("qvalue" %in% rownames(installed.packages())){
  library(qvalue)
} else {BiocManager::install("qvalue"); library(qvalue)}


CanCor <- function(x, y, xpfix=NULL, ypfix=NULL) 
 {
  if(is.null(xpfix)) { xpfix <- "x" }
  if(is.null(ypfix)) { ypfix <- "y" }
  tr <- function(r, n) r * sqrt((n-2)/(1-r^2))
  x <- scale(x)
  y <- scale(y)
  res <- cancor(x, y)
  ccres <- NULL
  ccres$cor$r <- res$cor
  ccres$cor$r[ccres$cor$r > 1] <- 1
  ccres$cor$r[ccres$cor$r < -1] <- -1
  names(ccres$cor$r) <- paste("CV", 1:length(ccres$cor$r), sep="")
  ccres[[xpfix]]$coef <- res$xcoef
  colnames(ccres[[xpfix]]$coef) <- paste(xpfix, "_CV", 1:ncol(ccres[[xpfix]]$coef), sep="")
  ccres[[ypfix]]$coef <- res$ycoef
  colnames(ccres[[ypfix]]$coef) <- paste(ypfix, "_CV", 1:ncol(ccres[[ypfix]]$coef), sep="")
  ccres[[xpfix]]$coef <- ccres[[xpfix]]$coef * sqrt(nrow(x) - 1)
  ccres[[ypfix]]$coef <- ccres[[ypfix]]$coef * sqrt(nrow(y) - 1)
  xy <- cbind(x[rownames(y),rownames(ccres[[xpfix]]$coef)], y[,rownames(ccres[[ypfix]]$coef)])
  n <- nrow(xy)
  xycor <- cor(xy)
  xn <- ncol(x)
  yn <- ncol(y)
  ccres$cor$rsq <- ccres$cor$r^2
  
  for(i in 1:length(ccres$cor$r)) {
   ccres$cor$chisq <- c(ccres$cor$chisq, -(n - 0.5 * (xn+yn+3)) * sum(log(1-ccres$cor$rsq[i:length(ccres$cor$r)])))
   ccres$cor$df <- c(ccres$cor$df, (xn-i+1)*(yn-i+1))
   ccres$cor$p <- c(ccres$cor$p, pchisq(ccres$cor$chisq[[i]], df=ccres$cor$df[[i]], lower.tail=FALSE))
  }
  if(all(class(try(qvalue(CCA_Pair$CV1p), silent = TRUE)) == "try-error")){
    ccres$cor$q <- qvalue(ccres$cor$p, pi0 = 1)$qvalues
  }else{
    ccres$cor$q <- qvalue(ccres$cor$p)$qvalues
  }
  ccres$cor <- t(data.frame(ccres$cor))
  
  ccres[[xpfix]]$score <- x[,rownames(ccres[[xpfix]]$coef)] %*% ccres[[xpfix]]$coef
  ccres[[ypfix]]$score <- y[,rownames(ccres[[ypfix]]$coef)] %*% ccres[[ypfix]]$coef
  
  ccres[[xpfix]]$loading$r <- xycor[rownames(ccres[[xpfix]]$coef),rownames(ccres[[xpfix]]$coef)] %*% ccres[[xpfix]]$coef
  ccres[[xpfix]]$loading$t <- tr(ccres[[xpfix]]$loading$r, n)
  ccres[[xpfix]]$loading$p <- pt(abs(ccres[[xpfix]]$loading$t), df=n-2, lower.tail=FALSE) * 2
  ccres[[xpfix]]$loading$p[ccres[[xpfix]]$loading$p > 1] <- 1

  ccres[[ypfix]]$loading$r <- xycor[rownames(ccres[[ypfix]]$coef),rownames(ccres[[ypfix]]$coef)] %*% ccres[[ypfix]]$coef
  ccres[[ypfix]]$loading$t <- tr(ccres[[ypfix]]$loading$r, n)
  ccres[[ypfix]]$loading$p <- pt(abs(ccres[[ypfix]]$loading$t), df=n-2, lower.tail=FALSE) * 2
  ccres[[ypfix]]$loading$p[ccres[[ypfix]]$loading$p > 1] <- 1

  ccres[[xpfix]]$importance$ContRatio <- colMeans(ccres[[xpfix]]$loading$r^2)
  ccres[[xpfix]]$importance$CumContRatio <- cumsum(ccres[[xpfix]]$importance$ContRatio)
  ccres[[xpfix]]$importance <- as.matrix(t(data.frame(ccres[[xpfix]]$importance)))
  
  ccres[[ypfix]]$importance$ContRatio <- colMeans(ccres[[ypfix]]$loading$r^2)
  ccres[[ypfix]]$importance$CumContRatio <- cumsum(ccres[[ypfix]]$importance$ContRatio)
  ccres[[ypfix]]$importance <- as.matrix(t(data.frame(ccres[[ypfix]]$importance)))
  
  ccres[[xpfix]]$crossloading$r <- xycor[rownames(ccres[[ypfix]]$coef),rownames(ccres[[xpfix]]$coef)] %*% ccres[[xpfix]]$coef
  ccres[[xpfix]]$crossloading$t <- tr(ccres[[xpfix]]$crossloading$r, n)
  ccres[[xpfix]]$crossloading$p <- pt(abs(ccres[[xpfix]]$crossloading$r), df=n-2, lower.tail=FALSE) * 2
  ccres[[xpfix]]$crossloading$p[ccres[[xpfix]]$crossloading$p > 1] <- 1
  
  ccres[[ypfix]]$crossloading$r <- xycor[rownames(ccres[[xpfix]]$coef),rownames(ccres[[ypfix]]$coef)] %*% ccres[[ypfix]]$coef
  ccres[[ypfix]]$crossloading$t <- tr(ccres[[ypfix]]$crossloading$r, n)
  ccres[[ypfix]]$crossloading$p <- pt(abs(ccres[[ypfix]]$crossloading$t), df=n-2, lower.tail=FALSE)
  ccres[[ypfix]]$crossloading$p[ccres[[ypfix]]$crossloading$p > 1] <- 1

  ccres[[xpfix]]$crossimportance$Redundancy <- colMeans(ccres[[xpfix]]$crossloading$r^2)
  ccres[[xpfix]]$crossimportance$CumRedundancy <- cumsum(ccres[[xpfix]]$crossimportance$Redundancy)
  ccres[[xpfix]]$crossimportance <- as.matrix(t(data.frame(ccres[[xpfix]]$crossimportance)))
  
  ccres[[ypfix]]$crossimportance$Redundancy <- colMeans(ccres[[ypfix]]$crossloading$r^2)
  ccres[[ypfix]]$crossimportance$CumRedundancy <- cumsum(ccres[[ypfix]]$crossimportance$Redundancy)
  ccres[[ypfix]]$crossimportance <- as.matrix(t(data.frame(ccres[[ypfix]]$crossimportance)))

  return(ccres)
}
