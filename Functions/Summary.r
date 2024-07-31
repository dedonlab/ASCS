############################################################
############################################################
# A modification of summary.gamlss function
# 
# Modified by: Shinsuke Ohnuki
# 
############################################################
############################################################
summary.gamlss <- function (object,
                           type = c("vcov", "qr"), robust = FALSE,
                           save = FALSE, 
                    hessian.fun = c("R", "PB"), 
                        verbose = FALSE,
                           ...) 
{
type <- match.arg(type)
if (type=="vcov")
 {
      covmat <- try(vcov(object, type="all", robust=robust,  hessian.fun = hessian.fun), silent = TRUE) 
             if (any(class(covmat)%in%"try-error"||any(is.na(covmat$se))))
                 { 
                   warning(paste("summary: vcov has failed, option qr is used instead\n"))
                   type <- "qr"
                 }
 }      
if (type=="vcov")
 {
        coef <- covmat$coef
          se <- covmat$se
      tvalue <- coef/se
      pvalue <-  2 * pt(-abs(tvalue), object$df.res)  
  coef.table <- cbind(coef, se, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef), c("Estimate" , "Std. Error" ,"t value","Pr(>|t|)"))  
      digits <- max(3, getOption("digits") - 3)
if(verbose) cat("*******************************************************************")
if(verbose) cat("\nFamily: ", deparse(object$family), "\n") 
if(verbose) cat("\nCall: ", deparse(object$call),  "\n", fill=TRUE)
if(verbose) cat("Fitting method:", deparse(object$method), "\n\n") 
    est.disp <- FALSE
#================ mu ESTIMATES ========================
    if ("mu"%in%object$parameters)   
    {
        if (object$mu.df != 0)   
        {          
                        pm <- object$mu.qr$rank 
                        p1 <- 1:pm
                        if(verbose) cat("-------------------------------------------------------------------\n")
                        if(verbose) cat("Mu link function: ", object$mu.link)
                        if(verbose) cat("\n")
                        if(verbose) cat("Mu Coefficients:")
            if (is.character(co <- object$contrasts)) 
              if(verbose) cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
                    paste, collapse = "="), "]")
            if(verbose) cat("\n")
            if(verbose) print.default(coef.table[p1,], digits = digits, print.gap = 2, quote = FALSE)
            if(verbose) cat("\n")
        }
        else
            if(object$mu.fix == TRUE) 
            {
              if(verbose) cat("-------------------------------------------------------------------\n")
              if(verbose) cat("Mu parameter is fixed \n")
            if (all(object$mu.fv == object$mu.fv[1]))
              if(verbose) cat("Mu = ", object$mu.fv[1], "\n")
            else
              if(verbose) cat("Mu is equal with the vector (", object$mu.fv[1], ",",object$mu.fv[2], ",",object$mu.fv[3], ",",object$mu.fv[4], ", ...) \n")
            }
    }
    if ("sigma"%in%object$parameters) 
    {
         if (object$sigma.df != 0)   
        {
                     ps <- object$sigma.qr$rank 
                     p1 <- (pm+1):(pm+ps)
        if(verbose) cat("-------------------------------------------------------------------\n")
        if(verbose) cat("Sigma link function: ", object$sigma.link)
        if(verbose) cat("\n")
        if(verbose) cat("Sigma Coefficients:")
        if(verbose) cat("\n")
        if(verbose) print.default(coef.table[p1,], digits = digits, print.gap = 2,  quote = FALSE)
        if(verbose) cat("\n")
        }
        else 
            if(object$sigma.fix == TRUE) 
            {
              if(verbose) cat("-------------------------------------------------------------------\n")
              if(verbose) cat("Sigma parameter is fixed")
              if(verbose) cat("\n")
            if (all(object$sigma.fv == object$sigma.fv[1]))
              if(verbose) cat("Sigma = ", object$sigma.fv[1], "\n")
            else
              if(verbose) cat("Sigma is equal with the vector (", object$sigma.fv[1], ",",object$sigma.fv[2], ",",object$sigma.fv[3], ",",object$sigma.fv[4], ", ...) \n")
            }
    }
    if ("nu"%in%object$parameters) 
    {
         if (object$nu.df != 0)   
        {
                   
                     pn <- object$nu.qr$rank 
                     p1 <- (pm+ps+1):(pm+ps+pn)
                     if(verbose) cat("-------------------------------------------------------------------\n")
                     if(verbose) cat("Nu link function: ", object$nu.link,"\n")
                     if(verbose) cat("Nu Coefficients:")
                     if(verbose) cat("\n")
                     if(verbose) print.default(coef.table[p1,], digits = digits, print.gap = 2,  quote = FALSE)
                     if(verbose) cat("\n")
        }
        else
            if(object$nu.fix == TRUE) {
              if(verbose) cat("-------------------------------------------------------------------\n")
              if(verbose) cat("Nu parameter is fixed")
              if(verbose) cat("\n")
            if (all(object$nu.fv == object$nu.fv[1]))
              if(verbose) cat("Nu = ", object$nu.fv[1], "\n")
            else
              if(verbose) cat("Nu is equal with the vector (", object$nu.fv[1], ",",object$nu.fv[2], ",",object$nu.fv[3], ",",object$nu.fv[4], ", ...) \n")
            }
    }
    if ("tau"%in%object$parameters) 
    {
         if (object$tau.df != 0)   
        {
                    
                     pt <- object$tau.qr$rank 
                     p1 <- (pm+ps+pn+1):(pm+ps+pn+pt)
                     if(verbose) cat("-------------------------------------------------------------------\n")
                     if(verbose) cat("Tau link function: ", object$tau.link,"\n")
                     if(verbose) cat("Tau Coefficients:")
                     if(verbose) cat("\n")
                     if(verbose) print.default(coef.table[p1,], digits = digits, print.gap = 2,  quote = FALSE)
                     if(verbose) cat("\n")
        }
        else
            if(object$tau.fix == TRUE) {
              if(verbose) cat("-------------------------------------------------------------------\n")
              if(verbose) cat("Tau parameter is fixed")
              if(verbose) cat("\n")
            if (all(object$tau.fv == object$tau.fv[1]))
              if(verbose) cat("Tau = ", object$tau.fv[1], "\n")
            else
              if(verbose) cat("Tau is equal with the vector (", object$tau.fv[1], ",",object$tau.fv[2], ",",object$tau.fv[3], ",",object$tau.fv[4], ", ...) \n")
            }
    }
if(verbose) cat("-------------------------------------------------------------------\n")
if(verbose) cat("No. of observations in the fit: ", object$noObs, "\n")
if(verbose)  cat("Degrees of Freedom for the fit: ", object$df.fit)
if(verbose) cat("\n")
if(verbose) cat("      Residual Deg. of Freedom: ", object$df.residual, "\n")
if(verbose) cat("                      at cycle: ", object$iter, "\n \n")
if(verbose) cat("Global Deviance:    ", object$G.deviance,
        "\n            AIC:    ",object$aic,
        "\n            SBC:    ",object$sbc, "\n")
if(verbose) cat("*******************************************************************")
if(verbose) cat("\n")
  } 
if (type=="qr")
  {
#-------------------------------------------------------------------------------
        estimatesgamlss<-function (object,Qr, p1, coef.p, est.disp , df.r, digits = max(3, getOption("digits") - 3),
                                               covmat.unscaled , ...)
          {
             dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
                                covmat <- covmat.unscaled
                                var.cf <- diag(covmat)
                                 s.err <- sqrt(var.cf)
                                tvalue <- coef.p/s.err
                                    dn <- c("Estimate", "Std. Error")
                          if (!est.disp) 
                             {
                                            pvalue <- 2 * pnorm(-abs(tvalue))
                                        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
                              dimnames(coef.table) <- list(names(coef.p), c(dn, "z value",
                                                    "Pr(>|z|)"))
                             }
                           else if (df.r > 0) 
                             {
                                           pvalue <- 2 * pt(-abs(tvalue), df.r)
                                       coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
                             dimnames(coef.table) <- list(names(coef.p), c(dn, "t value","Pr(>|t|)"))
                             }
                           else 
                           {
                                          coef.table <- cbind(coef.p, Inf)
                                dimnames(coef.table) <- list(names(coef.p), dn)
                           }
                          return(coef.table)
          }
      dispersion <- NULL
      digits <- max(3, getOption("digits") - 3)
if(verbose) cat("*******************************************************************")
if(verbose) cat("\nFamily: ", deparse(object$family), "\n") 
if(verbose) cat("\nCall: ", deparse(object$call),  "\n", fill=TRUE)
if(verbose) cat("Fitting method:", deparse(object$method), "\n\n") 
    est.disp <- FALSE
        df.r <- object$noObs - object$mu.df
#================ mu ESTIMATES ========================
    if ("mu"%in%object$parameters)   
    {
        if (object$mu.df != 0)   
        {
            Qr <- object$mu.qr 
            df.r <- object$noObs - object$mu.df
            if (is.null(dispersion)) 
                dispersion <- if (any(object$family == c("PO", 
                                    "BI", "EX", "P1"))) 
                                    1
                            else if (df.r > 0) {
                                                est.disp <- TRUE
                                                if (any(object$weights == 0)) 
                                                        warning(paste("observations with zero weight", 
                                                        "not used for calculating dispersion"))
                                                }
            else Inf   
            #---------------------------------------------------end of taken out 
                        p <- object$mu.df
                        p1 <- 1:(p-object$mu.nl.df)
            covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
            mu.coef.table <- estimatesgamlss(object=object,Qr=object$mu.qr, p1=p1, 
                                        coef.p=object$mu.coefficients[Qr$pivot[p1]], 
                                        est.disp =est.disp, df.r=df.r,
                                        covmat.unscaled =covmat.unscaled )
            if(verbose) cat("-------------------------------------------------------------------\n")
            if(verbose) cat("Mu link function: ", object$mu.link)
            if(verbose) cat("\n")
            if(verbose) cat("Mu Coefficients:")
            if (is.character(co <- object$contrasts)) 
              if(verbose) cat("  [contrasts: ", apply(cbind(names(co), co), 1, 
                    paste, collapse = "="), "]")
            if(verbose) cat("\n")
            if(verbose) print.default(mu.coef.table, digits = digits, print.gap = 2, quote = FALSE)
            if(verbose) cat("\n")
        }
        else
            if(object$mu.fix == TRUE) 
            {
              if(verbose) cat("-------------------------------------------------------------------\n")
              if(verbose) cat("Mu parameter is fixed")
              if(verbose) cat("\n")
            if (all(object$mu.fv == object$mu.fv[1]))
              if(verbose) cat("Mu = ", object$mu.fv[1], "\n")
            else
              if(verbose) cat("Mu is equal with the vector (", object$mu.fv[1], ",",object$mu.fv[2], ",",object$mu.fv[3], ",",object$mu.fv[4], ", ...) \n")
            }
        coef.table <- mu.coef.table
    }
    else
    {
        if (df.r > 0) {
                        est.disp <- TRUE
                        if (any(object$weights == 0)) 
                        warning(paste("observations with zero weight", 
                        "not used for calculating dispersion"))
                        }
    }
    if ("sigma"%in%object$parameters) 
    {
         if (object$sigma.df != 0)   
        {
                      Qr <- object$sigma.qr 
                    df.r <- object$noObs - object$sigma.df
                       p <- object$sigma.df 
                      p1 <- 1:(p-object$sigma.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        sigma.coef.table<-estimatesgamlss(object=object,Qr=object$sigma.qr, p1=p1, 
                                       coef.p=object$sigma.coefficients[Qr$pivot[p1]], 
                                       est.disp =est.disp, df.r=df.r,
                                       covmat.unscaled =covmat.unscaled )
        if(verbose) cat("-------------------------------------------------------------------\n")
        if(verbose) cat("Sigma link function: ", object$sigma.link)
        if(verbose) cat("\n")
        if(verbose) cat("Sigma Coefficients:")
        if(verbose) cat("\n")
        if(verbose) print.default(sigma.coef.table, digits = digits, print.gap = 2,  quote = FALSE)
        if(verbose) cat("\n")
        }
        else 
            if(object$sigma.fix == TRUE) {
              if(verbose) cat("-------------------------------------------------------------------\n")
              if(verbose) cat("Sigma parameter is fixed")
              if(verbose) cat("\n")
            if (all(object$sigma.fv == object$sigma.fv[1]))
              if(verbose) cat("Sigma = ", object$sigma.fv[1], "\n")
            else
              if(verbose) cat("Sigma is equal with the vector (", object$sigma.fv[1], ",",object$sigma.fv[2], ",",object$sigma.fv[3], ",",object$sigma.fv[4], ", ...) \n")
            } 
        if(!is.null(object$mu.qr))
          coef.table <- mu.coef.table
        if(!is.null(object$sigma.qr)) 
          coef.table <- rbind(coef.table, sigma.coef.table)
    }
    if ("nu"%in%object$parameters) 
    {
         if (object$nu.df != 0)   
        {
                     Qr <- object$nu.qr 
                   df.r <- object$noObs - object$nu.df
                      p <- object$nu.df 
                     p1 <- 1:(p-object$nu.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
          nu.coef.table <- estimatesgamlss(object=object,Qr=object$nu.qr, p1=p1, 
                                       coef.p=object$nu.coefficients[Qr$pivot[p1]], 
                                       est.disp =est.disp, df.r=df.r,
                                       covmat.unscaled =covmat.unscaled )
        if(verbose) cat("-------------------------------------------------------------------\n")
        if(verbose) cat("Nu link function: ", object$nu.link,"\n")
        if(verbose) cat("Nu Coefficients:")
        if(verbose) cat("\n")
        if(verbose) print.default(nu.coef.table, digits = digits, print.gap = 2,  quote = FALSE)
        if(verbose) cat("\n")
        }
        else
            if(object$nu.fix == TRUE) {
              if(verbose) cat("-------------------------------------------------------------------\n")
              if(verbose) cat("Nu parameter is fixed")
              if(verbose) cat("\n")
            if (all(object$nu.fv == object$nu.fv[1]))
              if(verbose) cat("Nu = ", object$nu.fv[1], "\n")
            else
              if(verbose) cat("Nu is equal with the vector (", object$nu.fv[1], ",",object$nu.fv[2], ",",object$nu.fv[3], ",",object$nu.fv[4], ", ...) \n")
            }
        if(!is.null(object$mu.qr))
          coef.table <- mu.coef.table
        if(!is.null(object$sigma.qr)) 
          coef.table <- rbind(coef.table, sigma.coef.table)
        if(!is.null(object$nu.qr)) 
          coef.table <- rbind(coef.table, nu.coef.table)
    }

    if ("tau"%in%object$parameters) 
   {
         if (object$tau.df != 0)   
        {
                     Qr <- object$tau.qr 
                   df.r <- object$noObs - object$tau.df
                      p <- object$tau.df 
                     p1 <- 1:(p-object$tau.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
         tau.coef.table <- estimatesgamlss(object=object,Qr=object$tau.qr, p1=p1, 
                                       coef.p=object$tau.coefficients[Qr$pivot[p1]], 
                                       est.disp =est.disp, df.r=df.r,
                                       covmat.unscaled =covmat.unscaled )
        if(verbose) cat("-------------------------------------------------------------------\n")
        if(verbose) cat("Tau link function: ", object$tau.link,"\n")
        if(verbose) cat("Tau Coefficients:")
        if(verbose) cat("\n")
        if(verbose) print.default(tau.coef.table, digits = digits, print.gap = 2,  quote = FALSE)
        if(verbose) cat("\n")
        }
        else
            if(object$tau.fix == TRUE) {
              if(verbose) cat("-------------------------------------------------------------------\n")
              if(verbose) cat("Tau parameter is fixed")
              if(verbose) cat("\n")
            if (all(object$tau.fv == object$tau.fv[1]))
              if(verbose) cat("Tau = ", object$tau.fv[1], "\n")
            else
              if(verbose) cat("Tau is equal with the vector (", object$tau.fv[1], ",",object$tau.fv[2], ",",object$tau.fv[3], ",",object$tau.fv[4], ", ...) \n")
            }
        if(!is.null(object$mu.qr))
          coef.table <- mu.coef.table
        if(!is.null(object$sigma.qr)) 
          coef.table <- rbind(coef.table, sigma.coef.table)
        if(!is.null(object$nu.qr)) 
          coef.table <- rbind(coef.table, nu.coef.table)
        if(!is.null(object$tau.qr)) 
          coef.table <- rbind(coef.table, tau.coef.table)
    }
if(verbose) cat("-------------------------------------------------------------------\n")
if(verbose) cat("No. of observations in the fit: ", object$noObs, "\n")
if(verbose) cat("Degrees of Freedom for the fit: ", object$df.fit)
if(verbose) cat("\n")
if(verbose) cat("      Residual Deg. of Freedom: ", object$df.residual, "\n")
if(verbose) cat("                      at cycle: ", object$iter, "\n \n")
if(verbose) cat("Global Deviance:    ", object$G.deviance,
        "\n            AIC:    ",object$aic,
        "\n            SBC:    ",object$sbc, "\n")
if(verbose) cat("*******************************************************************")
if(verbose) cat("\n")
    
  }
 if ( save == TRUE)
    {
      out <- as.list(environment())
       return(out)
  }
if(verbose) {
  invisible(coef.table)
} else {
  return(coef.table)
}
}