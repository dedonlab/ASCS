############################################################
############################################################
# Capturing combinatorial optimization
# 
# Author: Farzan Ghanegolmohammadi
# 
############################################################
############################################################
# Description:
## Detect influential factor(s) for differentiating functional clusters using logistic regression for
## canonical variables (CVs) as explanatory variables. 
############################################################
############################################################
setwd("PATH")
# Prerequisite
if("brglm" %in% rownames(installed.packages())){
  library(brglm)
} else {install.packages("brglm"); library(brglm)}

if("lmtest" %in% rownames(installed.packages())){
  library(lmtest)
} else {install.packages("lmtest"); library(lmtest)}

source("Functions/brstep.r")
############################################################
## Loading CCA results
load("Data/CCA/CCA.rdata")
## Loading HCA results
load("Data/GO/HCA.rdata")
## Significant level
FDR <- .05
############################################################
# Finding significant canonical variables (CV)
Sig_CVs <- colnames(CCA$cor)[CCA$cor["q",] < FDR]
############################################################
# Pairwise logistic regression
LogisticFit_Pairwise <- list()
Required_CV <- matrix(NA, nrow = length(HCA$CutTree_3), ncol = length(HCA$CutTree_3),
                      dimnames = list(names(HCA$CutTree_3),names(HCA$CutTree_3)))
Required_CV_LRT_P <- matrix(NA, nrow = length(HCA$CutTree_3), ncol = length(HCA$CutTree_3),
                            dimnames = list(names(HCA$CutTree_3),names(HCA$CutTree_3)))

pb <- txtProgressBar(min = 1, max = length(HCA$CutTree_3), style = 3)

for (i in names(HCA$CutTree_3)) {
  for (j in names(HCA$CutTree_3)[1:which(names(HCA$CutTree_3) == i)]) {
    if(j == i) next()
    ### Finding members of clusters in CV scores
    temp <- CCA$Codon$score[rownames(CCA$Codon$score) %in% c(HCA$CutTree_3[[i]], HCA$CutTree_3[[j]]), paste0("Codon_",Sig_CVs)]
    temp <- data.frame(Class = ifelse(rownames(temp) %in% HCA$CutTree_3[[i]], 1, 0), temp)
    temp <- temp[order(temp$Class),]

    # Logistic regression
    LogisticFit_Pairwise[[i]]$brglm <- brglm(Class ~ 1, data = temp, family =  binomial(logit), method = "brglm.fit",
                                           control.brglm = brglm.control(br.maxit = 2000L))
    ## Finding the best combination of the PC scores using step() function
    LogisticFit_Pairwise[[i]]$Step <- brstep(LogisticFit_Pairwise[[i]]$brglm, direction = "both", trace = 0,
                                           scope = as.formula(paste("~", paste(colnames(temp)[-1], collapse="+"), sep="")))
    if(length(coef(LogisticFit_Pairwise[[i]]$Step)) == 1){
      Required_CV[i,j] <- "Intercept"
    }else{
      reqCV <- names(coef(LogisticFit_Pairwise[[i]]$Step)[-1])
      reqCV <- strsplit(reqCV,"_")
      reqCV <- sapply(1:length(reqCV),function(x) reqCV[[x]][2])
      Required_CV[i,j] <- paste(reqCV, collapse = "|")
      rm(reqCV)
    }
    # Likelihood ratio test between the null and final models
    LogisticFit_Pairwise[[i]]$LRtest <- lrtest(LogisticFit_Pairwise[[i]]$brglm, LogisticFit_Pairwise[[i]]$Step)
    Required_CV_LRT_P[i,j] <- LogisticFit_Pairwise[[i]]$LRtest[2,"Pr(>Chisq)"]

    ### Estimating linear predictors
    LogisticFit_Pairwise[[i]]$Predictors <- predict(LogisticFit_Pairwise[[i]]$Step, newdata = temp, type = "link")

  }
  setTxtProgressBar(pb,which(names(HCA$CutTree_3) == i))
}
rm(i,j, pb)

save(LogisticFit_Pairwise, file = "Data/CCA/LogRegr_Pairwise.rdata")

write.csv(Required_CV, "Data/CCA/Pairwise_CV.csv")
write.csv(Required_CV_LRT_P, "Data/CCA/Pairwise_CV_LRT_P.csv.csv")

## q values
Required_CV_LRT_q <- qvalue::qvalue(Required_CV_LRT_P)$qvalues
write.csv(Required_CV_LRT_q, "Data/CCA/Pairwise_CV_LRT_q.csv.csv")