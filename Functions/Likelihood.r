############################################################
############################################################
# An additional function to gamlss package
# 
# Modefied by: Shinsuke Ohnuki
# 
############################################################
############################################################
# Description
## pdf function will switch by binomial family (BI, BB, ZIBI, ZIBB) to have bd argument.
############################################################
############################################################
gen.likelihood <- function(object)
{
  if (!is.gamlss(object)) stop("needs a gamlss object")
  fam <-  if(is.null(object$call$family)) as.gamlss.family(NO) 
  else as.gamlss.family(object$call$family)
  fname <- object$family[1]
  dfun <- paste("d",fname,sep="")
  pdf <- eval(parse(text=dfun))
  nopar <- length(object$par)
  y  <- object$y
  w  <- object$weights 
  X  <- list()
  links <- list()
  coefs <- list()
  Smo <- list() 
  offSet <- list()
  if (any(grepl("data", names(object$call)))) 
  {
    exitData <- TRUE  
    DaTa <- get(as.character(object$call["data"]))  
  }
  else
    exitData <- FALSE  
  if(fname == "BB" | fname == "BI" | fname == "ZIBI" | fname == "ZIBB") {
    isB <- TRUE
    y <- DaTa$y
    N <- DaTa$N
  }
  else {
    isB <- FALSE
  }
  for (i in object$par)
  {
    coefs[[i]] <- eval(parse(text=paste(paste("object$", i, sep=""), ".coefficients", sep="")))
    notNAcoef <- !is.na(coefs[[1]])
    X[[i]] <- eval(parse(text=paste(paste("object$", i, sep=""), ".x", sep="")))  
    if (any(is.na(coefs[[1]])))
    {
      coefs[[i]] <- coefs[[i]][notNAcoef]
      X[[i]] <- as.matrix(X[[i]][,notNAcoef])
    }
    links[i] <- eval(parse(text=paste(paste("object$", i, sep=""), ".link", sep="")))
    Smo[[i]] <- eval(parse(text=paste(paste("object$", i, sep=""), ".s", sep="")))
    offSet[[i]] <- eval(parse(text=paste(paste("object$", i, sep=""), ".offset", sep=""))) 
  } 
  switch(nopar,
{ #  1 parameter
  lik.fun <- function(par) 
  {
    lmu <-  length(coefs[["mu"]])
    if (length(par)!=lmu) stop("par is not the right length")
    eta.mu <- if (is.null(Smo[["mu"]])) X[["mu"]] %*% par[1:lmu]   + offSet[["mu"]]  
    else X[["mu"]] %*% par[1:lmu] + rowSums(Smo[["mu"]]) + offSet[["mu"]] 
    mu <- fam$mu.linkinv(eta.mu) 
    -sum(w*pdf(y, mu=mu, log=TRUE))
  }
  thebetas <-  coefs[["mu"]] 
  formals(lik.fun) <- alist(par=thebetas)
},
{ # 2 parameter
  lik.fun <- function(par) 
  { 
    lmu <- length(coefs[["mu"]])
    lsigma <- length(coefs[["sigma"]])
    tl <- lmu + lsigma	
    if (length(par)!=tl) stop("par is not the right length")	
    eta.mu <- if (is.null(Smo[["mu"]])) X[["mu"]] %*% par[1:lmu]  + offSet[["mu"]] 
    else X[["mu"]] %*% par[1:lmu]+ rowSums(Smo[["mu"]]) + offSet[["mu"]] 
    mu <- fam$mu.linkinv(eta.mu) 
    eta.sigma <- if (is.null(Smo[["sigma"]])) X[["sigma"]] %*% par[(lmu+1):(lmu+lsigma)] + offSet[["sigma"]] 
    else X[["sigma"]] %*% par[(lmu+1):(lmu+lsigma)] + rowSums(Smo[["sigma"]])+offSet[["sigma"]]
    sigma <- fam$sigma.linkinv(eta.sigma)
    if(isB) {
      -sum(w*pdf(y, mu=mu, sigma=sigma, bd=N, log=TRUE))
    } else {
      -sum(w*pdf(y, mu=mu, sigma=sigma, log=TRUE))
    }
  }
  thebetas <-  c(coefs[["mu"]], coefs[["sigma"]]) 
  formals(lik.fun) <- alist(par=thebetas)
},         
{ # 3 parameter 
  lik.fun <- function(par) 
  {
    lmu <-	 length(coefs[["mu"]])
    lsigma <-  length(coefs[["sigma"]])
    lnu <-  length(coefs[["nu"]])
    tl <-  lmu + lsigma + lnu	
    if (length(par)!=tl) stop("par is not the right length")	
    eta.mu <- if (is.null(Smo[["mu"]])) X[["mu"]] %*% par[1:lmu] + offSet[["mu"]] 
    else X[["mu"]] %*% par[1:lmu]+ rowSums(Smo[["mu"]]) + offSet[["mu"]] 
    mu <- fam$mu.linkinv(eta.mu) 
    eta.sigma <- if (is.null(Smo[["sigma"]])) X[["sigma"]] %*% par[(lmu+1):(lmu+lsigma)] +offSet[["sigma"]] 
    else X[["sigma"]] %*% par[(lmu+1):(lmu+lsigma)] + rowSums(Smo[["sigma"]])+ offSet[["sigma"]]       
    sigma <- fam$sigma.linkinv(eta.sigma)
    eta.nu <- if (is.null(Smo[["nu"]])) X[["nu"]] %*% par[(lmu+lsigma+1):(lmu+lsigma+lnu)] + offSet[["nu"]] 
    else X[["nu"]] %*% par[(lmu+lsigma+1):(lmu+lsigma+lnu)] + rowSums(Smo[["nu"]]) + offSet[["nu"]]   
    nu <- fam$nu.linkinv(eta.nu)     
    -sum(w*pdf(y, mu=mu, sigma=sigma, nu=nu, log=TRUE)) 	
  }
  thebetas <-  c(coefs[["mu"]], coefs[["sigma"]], coefs[["nu"]]) 
  formals(lik.fun) <- alist(par=thebetas)
},
{ # 4 parameter
  lik.fun <- function(par) 
  {
    lmu <- length(coefs[["mu"]])
    lsigma <- length(coefs[["sigma"]])
    lnu <-  length(coefs[["nu"]])
    ltau <-  length(coefs[["tau"]])
    tl <- lmu + lsigma + lnu +ltau	
    if (length(par)!=tl) stop("par is not the right length")	
    eta.mu <- if (is.null(Smo[["mu"]])) X[["mu"]] %*% par[1:lmu]  + offSet[["mu"]] 
    else X[["mu"]] %*% par[1:lmu]+ rowSums(Smo[["mu"]]) + offSet[["mu"]] 
    mu <- fam$mu.linkinv(eta.mu) 
    eta.sigma <- if (is.null(Smo[["sigma"]])) X[["sigma"]] %*% par[(lmu+1):(lmu+lsigma)] + offSet[["sigma"]]  
    else X[["sigma"]] %*% par[(lmu+1):(lmu+lsigma)] + rowSums(Smo[["sigma"]]) + offSet[["sigma"]] 
    sigma <- fam$sigma.linkinv(eta.sigma)
    eta.nu <- if (is.null(Smo[["nu"]])) X[["nu"]] %*% par[(lmu+lsigma+1):(lmu+lsigma+lnu)]+ offSet[["nu"]] 
    else X[["nu"]] %*% par[(lmu+lsigma+1):(lmu+lsigma+lnu)] + rowSums(Smo[["nu"]])+ offSet[["nu"]]     
    nu <- fam$nu.linkinv(eta.nu) 
    eta.tau <- if (is.null(Smo[["tau"]]))  X[["tau"]] %*% par[(lmu+lsigma+lnu+1):(lmu+lsigma+lnu+ltau)] + offSet[["tau"]] 
    else  X[["tau"]] %*% par[(lmu+lsigma+lnu+1):(lmu+lsigma+lnu+ltau)]+ rowSums(Smo[["tau"]]) + offSet[["tau"]] 
    tau <- fam$tau.linkinv(eta.tau)
    -sum(w*pdf(y, mu=mu, sigma=sigma, nu=nu, tau=tau, log=TRUE)) 	
  }
  thebetas <-  c(coefs[["mu"]], coefs[["sigma"]], coefs[["nu"]], coefs[["tau"]]) 
  formals(lik.fun) <- alist(par=thebetas)
})
lik.fun
}