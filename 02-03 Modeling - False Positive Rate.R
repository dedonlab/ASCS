############################################################
############################################################
# Estimating false positive rate
# 
# Author: Farzan Ghanegolmohammadi
#
############################################################
############################################################
# Description:
## Implementation of parametric bootstrapping.
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

if("foreach" %in% rownames(installed.packages())){
  library(foreach)
} else {install.packages("foreach"); library(foreach)}
## For Windows
if("doSNOW" %in% rownames(installed.packages())){
  library(doSNOW)
} else {install.packages("doSNOW"); library(doSNOW)}

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
## Loading fitted model (ICF ~ GC3)
load("Data/GLM/BIfit.rdata")
# Reading GC3 values
GC3 <- read.csv("Data/GC3.csv", header = T, row.names = 1)
## Number of iterations
itr <- 2000
############################################################
# Parametric bootstrapping
cl <- parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(cl)


foreach(j = 1:itr, .packages = c("gamlss","qvalue"), .inorder = TRUE) %dopar% {
          gcont <- gamlss.control(n.cyc = 200L, trace = FALSE)

          ## A matrix to save p-values
          BIPValues_Rand <-  matrix(data = NA ,nrow = nrow(CodonCount_N), ncol = ncol(CodonCount_N),
                                    dimnames = list(rownames(CodonCount_N), colnames(CodonCount_N)))

          for (i in C59) {
            ## 1. Getting mu from fitted model
            mu_hat <- fitted(BIfit[[i]], what = "mu")

            ## 2. Simulating codon counts under binomial model
            n_sim <- rBI(n = nrow(CodonCount_N), bd = CodonCount_N[,i], mu = mu_hat)

            ## 3. Fitting model to simulated data
            temp <- data.frame(n = n_sim, N = CodonCount_N[,i], x = GC3[rownames(CodonCount_n),"GC3"],
                               row.names = rownames(CodonCount_N))
            fit_sim <- gamlss(cbind(n, N-n) ~ x, data = temp, family = BI, control = gcont)

            ## 4. Getting gene-specific mu from the simulated model
            mu_sim <- fitted(fit_sim, what = "mu")
            names(mu_sim) <- rownames(CodonCount_N)

            ## A matrix to save P values of lower- and upper-tail
            res <- matrix(NA, nrow = nrow(temp), ncol = 2, dimnames = list(rownames(temp), c("LowerTail","UpperTail")))

            ## 5. Computing p-values (2-sided)
            for (k in rownames(temp)) {
              if(temp[k,"n"] == 0 & temp[k,"N"] == 0){ # When a given AA does not exist in the protein sequence (n == 0 & N == 0)
                res[k,"LowerTail"] <- NA
                res[k,"UpperTail"] <- NA
                next()
              }
              ### Lower tail
              res[k,"LowerTail"] <- pBI(q = temp[k,"n"], bd = temp[k,"N"], mu = mu_sim[k], lower.tail = T)
              ### Upper tail
              ## Correcting for n == 0 or  n == N
              if(temp[k,"n"] == 0 | temp[k,"n"] == temp[k,"N"]){
                ### Because of the coding bugs, when n == 0 or  n == N and lower.tail = F; then, pBI() + dBI() should be used
                res[k,"UpperTail"] <- pBI(q = temp[k,"n"], bd = temp[k,"N"], mu = mu_sim[k], lower.tail = F) +
                                      dBI(x = temp[k,"n"], bd = temp[k,"N"], mu = mu_sim[k], log = F)
              }else{ ## when n!= 0 or n != N
                res[k,"UpperTail"] <- pBI(q = temp[k,"n"], bd = temp[k,"N"], mu = mu_sim[k], lower.tail = F)
              }
            }
            ### Two-sided test
            P2Sided <- 2 * (apply(res, MARGIN = 1, min))
            BIPValues_Rand[names(P2Sided),i] <- P2Sided
            rm(temp, P2Sided, k, res)
          }
          ## Saving the results in each iteration
          save(BIPValues_Rand, file = paste0("Data/GLM/BIPValues_Rand(",sprintf("Rdm%04i",j),").rdata"))
        }

stopCluster(cl)
###################################
# Estimating false positive rate
SigORFs_Rand <- list()
pb <- txtProgressBar(min = 1, max = itr, style = 3)
for (j in 1:itr) {
  load(paste0("Data/GLM/BIPValues_Rand(",sprintf("Rdm%04i",j),").rdata"))

  BIPValues_Rand[BIPValues_Rand > 1] <- 1
  BIqValues_GC3_Rand <- qvalue::qvalue(BIPValues_Rand[,C59])$qvalues
  BIqValues_GC3_Rand[is.na(BIqValues_GC3_Rand)] <- 1
  SigORFs_Rand$FDR_10[[as.character(j)]] <- rownames(BIqValues_GC3_Rand)[rowSums(BIqValues_GC3_Rand < 0.1, na.rm = T) > 0]
  SigORFs_Rand$FDR_05[[as.character(j)]] <- rownames(BIqValues_GC3_Rand)[rowSums(BIqValues_GC3_Rand < 0.05, na.rm = T) > 0]
  SigORFs_Rand$FDR_01[[as.character(j)]] <- rownames(BIqValues_GC3_Rand)[rowSums(BIqValues_GC3_Rand < 0.01, na.rm = T) > 0]
  rm(BIPValues_Rand);gc()
  setTxtProgressBar(pb, which(1:itr == j))
}
rm(j)
save(SigORFs_Rand, file = "Data/GLM/SigORFs_Rand.rdata")
