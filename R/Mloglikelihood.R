################################################################################
#  Minus the log-likelihood                                                    #
################################################################################
#                                                                              #
#  This function computes minus the logarithm of the likelihood function       #
#                                                                              #
#  Its parameters are                                                          #
#   - p      : the parameters vector, in the form                              #
#              c( frailty distribution parameter(s),                           #
#                 baseline hazard parameter(s),                                #
#                 regression parameter(s) )                                    #
#   - obs    : the observed data, in the form                                  #
#              list( time   = event/censoring times,                           #
#                   [trunc  = left truncation times, ]                         #
#                    event   = event indicators,                               #
#                    x       = covariate data.frame, intercept included        #
#                    cluster = cluster ID vector,                              #
#                    ncl     = number of clusters,                             #
#                    di      = vector giving the numbers of events per cluster #
#   - dist   : the baseline hazard distribution name                           #
#   - frailty: the frailty distribution name                                   #
#   - correct  : (only for possta) the correction to use in case of many       #
#                events per cluster to get finite likelihood values.           #
#                When correct!=0 the likelihood is divided by                  #
#                10^(#clusters * correct) for computation,                     #
#                but the value of the log-likelihood in the output             #
#                is the re-adjusted value.                                     #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
#   Last modification on: February 1, 2012                                     #
################################################################################

Mloglikelihood <- function(p,
                           obs,
                           dist,
                           frailty,
                           correct) { 
  # ---- Assign the number of frailty parameters 'nFpar' ---------------------#
  # ---- and compute Omega for the Positive Stable frailty -------------------#
  
  if (frailty == "none") 
    nFpar <- 0
  else if (frailty %in% c("gamma","ingau")) {
    theta <- exp(p[1])
    nFpar <- 1 
  }
  else if (frailty == "possta") {
    nu <- exp(-exp(p[1]))
    nFpar <- 1
    D <- max(obs$di)
    Omega <- Omega(D, correct=correct, nu=nu)
  }
  
    
  # ---- Baseline hazard -----------------------------------------------------#
  
  # baseline parameters
  if (dist == "weibull") { 
    pars <- c(rho    = exp(p[nFpar+1]),
              lambda = exp(p[nFpar+2]))
    beta <- p[-(1:(nFpar+2))]
  } else
  if (dist == "exponential") {
    pars <- c(lambda = exp(p[nFpar+1]))
    beta <- p[-(1:(nFpar+1))]
  } else
  if (dist == "gompertz") {
    pars <- c(gamma  = exp(p[nFpar+1]),
              lambda = exp(p[nFpar+2]))    
    beta <- p[-(1:(nFpar+2))]
  } else
  if (dist == "lognormal") {
    pars <- c(mu    = p[nFpar+1],
              sigma = exp(p[nFpar+2]))    
    beta <- p[-(1:(nFpar+2))]
  } else
  if (dist == "loglogistic") {
    pars <- c(alpha = p[nFpar+1],
              kappa = exp(p[nFpar+2]))    
    beta <- p[-(1:(nFpar+2))]
  }
  
  # baseline: from string to the associated function
  dist <- eval(parse(text=dist))

    
  # ---- Cumulative Hazard by cluster ----------------------------------------#
  
  cumhaz <- NULL
  if (frailty == "none") { ### NO FRAILTY
    cumhaz <- sum(dist(pars, obs$time, what="H") * 
      exp(as.matrix(obs$x) %*% c(0, beta)))
   
    # Possible truncation
    if (!is.null(obs$trunc)) 
      cumhaz <- cumhaz - 
        sum(dist(pars, obs$trunc, what="H") * 
        exp(as.matrix(obs$x) %*% c(0, beta)))
  } else { ### FRAILTY
  cumhaz <- aggregate(
    dist(pars, obs$time, what="H") *
    exp(as.matrix(obs$x) %*% c(0, beta)), 
    by=list(obs$cluster), FUN=sum)[, 2]
    
    # Possible truncation
    if (!is.null(obs$trunc)) 
      cumhaz <- cumhaz - aggregate(
        dist(pars, obs$trunc, what="H") *
        exp(as.matrix(obs$x) %*% c(0, beta)), 
        by=list(obs$cluster), FUN=sum)[, 2]
  } 
    
  # ---- log-hazard by cluster -----------------------------------------------#
  
  loghaz <- NULL
  if (frailty == "none") {
    loghaz <- sum(obs$event * (dist(pars, obs$time, what="lh") + 
      as.matrix(obs$x) %*% c(0, beta)))
  } else {
    loghaz <- aggregate(obs$event * (dist(pars, obs$time, what="lh") + 
      as.matrix(obs$x) %*% c(0, beta)),
                        by=list(obs$cluster), FUN=sum)[, 2]
  }

    
  # ---- log[ (-1)^di L^(di)(cumhaz) ]---------------------------------------#
    
  logSurv <- NULL
  if (frailty=="none")
    logSurv <- mapply(fr.none, s=cumhaz, what="logLT") else
  if (frailty=="gamma")
    logSurv <- mapply(fr.gamma, 
                      k=obs$di, s=cumhaz, theta=rep(theta, obs$ncl), 
                      what="logLT") 
  else if (frailty=="ingau")
    logSurv <- mapply(fr.ingau, 
                      k=obs$di, s=cumhaz, theta=rep(theta, obs$ncl), 
                      what="logLT") 
  else if (frailty=="possta"){
    logSurv <- sapply(1:obs$ncl, 
                      function(x) fr.possta(k=obs$di[x], s=cumhaz[x], 
                                            nu=nu, Omega=Omega, 
                                            what="logLT",
                                            correct=correct))
  }

    
  # ---- Minus the log likelihood --------------------------------------------#
  
  Mloglik <- -sum(loghaz + logSurv)
  attr(Mloglik, "cumhaz") <- cumhaz
  return(Mloglik)
}

