############################################################
############################################################
# Codon extraction
# 
# Author: Farzan Ghanegolmohammadi
#
############################################################
############################################################
# Description:
## The ORFs of Saccharomyces cerevisiae S288C strain (i.e., reference strain)
## was downloaded from The Saccharomyces Genome Database (SGD; https://www.yeastgenome.org/).
## Blocked reading frames and transposable elements were deleted from nucleus genome.
############################################################
############################################################
setwd("PATH/Data/")

# Prerequisite
if("seqinr" %in% rownames(installed.packages())){
  library(seqinr)
} else {install.packages("seqinr"); library(seqinr)}
############################################################
# Codon information
CodonInfo <- read.csv("Data/CodonInfo.csv", header = TRUE, row.names = 1)
## 59 codons
C59 <- CodonInfo$Codon[!CodonInfo$Codon %in% c("ATG", "TGG", "TAA","TAG","TGA")]
# Reading Coding Region Sequences (CDS) of yeast (S. cerevisiae strain S288C)
CDS <- read.fasta(file = "Data/SGD/orf_coding_5797.fasta", seqtype = "DNA",
                  forceDNAtolower = FALSE, set.attributes = FALSE)
############################################################
# Count data
CodonCount_n <- matrix(NA, nrow = length(CDS), ncol = nrow(CodonInfo),
                       dimnames = list(names(CDS), CodonInfo$Codon))
pb <- txtProgressBar(min = 0, max = length(CDS), style = 3)

for (i in names(CDS)) {
  ll <- seq(1, length(CDS[[i]]),3)
  cc <- sapply(ll, function(x) paste(CDS[[i]][x:(x+2)],collapse = ""))
  if(length(cc) * 3 != length(CDS[[i]])) stop()
  res <- table(cc) 
  CodonCount_n[i, names(res)] <- res
  rm(ll,cc,res)
  setTxtProgressBar(pb,which(names(CDS) == i))
}
rm(i,pb)

# Replace NAs with 0
CodonCount_n[is.na(CodonCount_n)] <- 0

write.csv(CodonCount_n, "Data/CodonCount_n.csv")
############################################################
# Sum of synonymous codons
CodonCount_N <- matrix(NA, nrow = length(CDS), ncol = nrow(CodonInfo),
                       dimnames = list(names(CDS), CodonInfo$Codon))

pb <- txtProgressBar(min = 0, max = length(CDS), style = 3)

for (i in names(CDS)) {
  for(j in unique(CodonInfo$AA)){
    cc <- CodonInfo[CodonInfo$AA %in% j,"Codon"]
    CodonCount_N[i, cc] <- sum(CodonCount_n[i, cc], na.rm = T)
    rm(cc)
  }
  setTxtProgressBar(pb,which(names(CDS) == i))
}
rm(i,j,pb)

write.csv(CodonCount_N, "Data/CodonCount.csv")
