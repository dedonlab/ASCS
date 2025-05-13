############################################################
############################################################
# GC3 estimation
# 
# Author: Farzan Ghanegolmohammadi
#
############################################################
############################################################
# Description:
## Estimating GC3 metric using formula "3" Sahoo et al. (2019; Gene; https://doi.org/10.1016/j.gene.2019.100012)
############################################################
############################################################
setwd("PATH/")

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
# Estimating GC3
GC3 <- rep(NA, length(CDS))
names(GC3) <- names(CDS)

for (i in names(CDS)) {
  ll <- seq(1, length(CDS[[i]]), 3)
  codon_seq <- sapply(ll, function(x) toupper(paste(CDS[[i]][x:(x+2)],collapse = "")))
  codon_seq <- codon_seq[codon_seq %in% C59]
  if (length(codon_seq) == 0) {
    GC3[i] <- NA
  } else {
    GC3[i] <- sum(grepl("[GC]$", codon_seq)) / length(codon_seq)
  }
}

write.csv(data.frame(GC3 = GC3), "GC3.csv", row.names = T)
