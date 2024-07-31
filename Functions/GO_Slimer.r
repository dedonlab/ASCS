############################################################
############################################################
# Gene Ontology (GO) slimmer
# 
# Author: Shinsuke Ohnuki
# 
############################################################
############################################################
goSlimer <- function(gomx)
 {
  temp <- gomx
  res <- NULL
  pb <- txtProgressBar(min = 0, max = ncol(gomx), style = 3)
  for(i in colnames(gomx)) {
   if(sum(colnames(temp) %in% i) > 0) {
    res[[i]] <- i
    if(!is.list(res)) res <- as.list(res)
    temp <- temp[,colnames(temp) != i,drop=FALSE]
    for(j in colnames(temp)) {
     #browser()
     if(sum(gomx[,i] == temp[,j]) == nrow(gomx)) res[[i]] <- c(res[[i]], j)
    }
    #browser()
    if(sum(!colnames(temp) %in% res[[i]]) > 0) {
     temp <- temp[,!colnames(temp) %in% res[[i]], drop=FALSE]
    }
   }
    setTxtProgressBar(pb, which(colnames(gomx) == i))
  }
  close(pb)
  return(res)
 }
