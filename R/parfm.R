################################################################################
#  Parametric frailty models fitting                                           #
################################################################################
#                                                                              #
#  This function is the core of the package,                                   #
#   performing estimation and returning results.                               #
#                                                                              #
#  It is also the front-end of the package,                                    #
#   being the only that will be visible to the user.                           #
#                                                                              #
#  Its parameters are                                                          #
#   - formula  : a formula object, with the response                           #
#                on the left of a ~ operator,                                  #
#                and the terms on the right.                                   #
#                The response must be a survival object                        #
#                as returned by the Surv() function;                           #
#   - cluster  : the name of the variable in data containing cluster IDs;      #
#   - data     : a data.frame containing all the variables;                    #
#   - inip     : initial values for the baseline hazard & reg parameters;      #
#   - iniFpar  : initial value(s) for the heterogeneity parameter(s);          #
#   - dist     : the baseline hazard;                                          #
#   - frailty  : the frailty distribution;                                     #
#   - method   : the optimisation method (See optim());                        #
#   - maxit    : the maximum number of iteration (See optim());                #
#   - Fparscale: the scaling value for all the frailty parameter(s) in optim() #
#                Optimisation is performed on Fpar/Fparscale                   #
#   - showtime : is the execution time displayed?                              #
#   - correct  : (only for possta) the correction to use in case of many       #
#                events per cluster to get finite likelihood values.           #
#                When correct!=0 the likelihood is divided by                  #
#                10^(#clusters * correct) for computation,                     #
#                but the value of the log-likelihood in the output             #
#                is the re-adjusted value.                                     #
#                                                                              #
#                                                                              #
#   Date: December 21, 2011                                                    #
#   Last modification on: February 15, 2012                                    #
################################################################################

parfm <- function(formula,
                  cluster=NULL,
                  data,
                  inip=NULL,
                  iniFpar=NULL,
                  dist="weibull",
                  frailty="none",
                  method="BFGS",
                  maxit=500,
                  Fparscale=1,
                  showtime=TRUE,
                  correct=0){

  #----- Check the baseline hazard and the frailty distribution ---------------#
  if (!(dist %in% 
    c("exponential", "weibull", "gompertz", "loglogistic", "lognormal")))
    stop("invalid baseline hazard")
  if (!(frailty %in% 
    c("none", "gamma", "ingau", "possta")))
    stop("invalid frailty distribution")
  if (frailty == "none" &&  !is.null(cluster))
    warning(paste("With frailty='none' the cluster variable '",
                  cluster, "' is not used!", sep=""))
  if (frailty == "none" &&  !is.null(iniFpar))
    warning("With frailty='none' the argument 'iniFpar' is not used!")
  
  #----- 'Correct' is useless except for frailty="possta" ---------------------#
  if (frailty == "possta") {  #Do not exaggerate when setting 'correct' !
    if (10^correct == Inf || 10^-correct == 0)
      stop("'correct' is too large!")
    if (10^correct == 0 || 10^-correct == Inf)
      stop("'correct' is too small!")
  }
  else if (correct != 0)
      warning(paste("'correct' has no effect when 'frailty = ", frailty, "'",
                    sep=""))
  
  #----- Data for Mloglikelihood() --------------------------------------------#
  obsdata <- NULL
    #time
  if (length(formula[[2]]) == 3) {          # --> without left truncation
    obsdata$time <- eval(parse(text=paste("data$", 
                                          formula[[2]][[2]], sep="")))
    obsdata$event <- eval(parse(text=paste("data$", 
                                          formula[[2]][[3]], sep=""))) 
  } else if (length(formula[[2]]) == 4) {   # --> with left truncation
    obsdata$trunc <- eval(parse(text=paste("data$", 
                                          formula[[2]][[2]], sep="")))
    
    obsdata$time <- eval(parse(text=paste("data$", 
                                          formula[[2]][[3]], sep="")))
    obsdata$event <- eval(parse(text=paste("data$", 
                                          formula[[2]][[4]], sep="")))
  }
  if (!all(levels(as.factor(obsdata$event)) %in% 0:1))
  stop(paste("The status indicator 'event' in the Surv object",
              "in the left-hand side of the formula object",
              "must be either 0 (no event) or 1 (event)."))

    #covariates (an intercept is automatically added)
  obsdata$x <- as.data.frame(model.matrix(formula, data=data))
    #cluster
  if (is.null(cluster)) {
    if (frailty != "none")
      stop(paste("if you specify a frailty distribution,\n",
                 "then you have to specify the cluster variable as well"))
  } else {
    if (! cluster %in% names(data))
      stop(paste("object '", cluster, "' not found", sep=""))
    obsdata$cluster <- eval(parse(text=paste("data$", cluster, sep="")))
      #number of clusters
    obsdata$ncl <- length(levels(as.factor(obsdata$cluster)))
      #number of events in each cluster
    obsdata$di <- aggregate(obsdata$event, 
                            by=list(obsdata$cluster), FUN=sum)[, 2]
  }
  
  #----- Dimensions -----------------------------------------------------------#
    #nFpar: number of heterogeneity parameters
  if (frailty == "none")
    nFpar <- 0
  else if (frailty %in% c("gamma", "ingau", "possta"))
    nFpar <- 1
  
    #nBpar: number of parameters in the baseline hazard
  if (dist == "exponential")
    nBpar <- 1
  else if (dist %in% c("weibull", "gompertz", "lognormal", "loglogistic"))
    nBpar <- 2
  
    #nRpar: number of regression parameters
  nRpar <- ncol(obsdata$x) - 1
  
  #----- Initial parameters ---------------------------------------------------#
  if (!is.null(inip)) {
      #if they are specified, then 
      #(1) check the dimension,
      #(2) check whether they lie in their parameter space, and
      #(3) reparametrise them so that they take on values on the whole real line
    
    if (length(inip) != nBpar + nRpar)
      stop(paste("number of initial parameters 'inip' must be", nBpar + nRpar))
    p.init <- inip
    if (dist %in% c("exponential", "weibull", "gompertz"))
      #1st initial par: log(lambda), log(rho), or log(gamma)
      if (p.init[1] <= 0)
        stop(paste("with that baseline, the 1st parameter has to be > 0"))
      p.init[1] <- log(p.init[1]) 
    if (dist %in% c("weibull", "gompertz", "lognormal", "loglogistic")) {
      #2nd initial par: log(lambda), log(lambda), log(sigma), or log(kappa)
      if (p.init[2] <= 0)
        stop(paste("with that baseline, the 2nd parameter has to be > 0"))
      p.init[2] <- log(p.init[2]) 
    }
  } else {
      #if they are not specified, then fit a parametric Cox's model
    require(eha)
    
    shape <- 0  #if zero or negative, the shape parameter is estimated
    d <- dist
    if (d == "exponential") {
      d <- "weibull"
      shape <- 1  #if positive, the shape parameter is fixed at that value
    }
    
    coxMod <- phreg(formula=formula, data=data,
                    dist=d, shape=shape, 
                    center=FALSE, control=list(maxiter=maxit))
    logshape <- as.numeric(coxMod$coef["log(shape)"])
    logscale <- as.numeric(coxMod$coef["log(scale)"])
    
    if (dist == "exponential") {
      p.init <- - logscale                        #log(lambda)
    } 
    else if (dist == "weibull") {
      p.init <- c(logshape,                       #log(rho)
                  - exp(logshape) * logscale)     #log(lambda)
    } 
    else if (dist == "gompertz") {
      p.init <- c(- logscale,                     #log(gamma)
                  logshape)                       #log(lambda)
    } 
    else if (dist == "lognormal") {
      p.init <- c(logscale,                       #mu
                  - logshape)                     #log(sigma)
    } 
    else if (dist == "loglogistic") {
      p.init <- c(- exp(logshape) * logscale,     #alpha
                  logshape)                       #log(kappa)
    }
    
    if (nRpar > 0)
      p.init <- c(p.init,
                  as.numeric(coxMod$coef[1:nRpar]))   
  }
  
  #--- frailty parameters initialisation ---#
  if (frailty == "none")
    pars <- NULL
  else if (frailty %in% c("gamma", "ingau")) {
    if (is.null(iniFpar))
      iniFpar <- 1
    else if (iniFpar <= 0)
      stop("initial heterogeneity parameter (theta) has to be > 0")
    pars <- log(iniFpar)
  }
  else if (frailty == "possta") {
    if (is.null(iniFpar))
      iniFpar <- 0.5
    else if (iniFpar <= 0 || iniFpar >= 1)
      stop("initial heterogeneity parameter (nu) must lie in (0, 1)")
    pars <- log(-log(iniFpar))
  }
  
  pars <- c(pars, p.init)
  
  #----- Minimise Mloglikelihood() --------------------------------------------#
  if ((frailty == "none") && is.null(inip)) {
    todo <- expression({res <- list(par=pars)})
    if (showtime)
      extime <- system.time(eval(todo))[1]
    else {
      eval(todo)
      extime <- NULL
    }
 
    it <- NULL
    lL <- coxMod$loglik[2]
  }
  else {
    todo <- expression({
      res <- optim(par=pars, fn=Mloglikelihood, method=method, 
                   obs=obsdata, dist=dist, frailty=frailty,
                   correct=correct,
                   hessian=TRUE, 
                   control=list(maxit=maxit, 
                                parscale=c(rep(Fparscale, nFpar),
                                           rep(1, nBpar+nRpar))))})
    if (showtime)
      extime <- system.time(eval(todo))[1]
    else {
      eval(todo)
      extime <- NULL
    }
 
    if(res$convergence > 0)
      warning("optimisation procedure did not converge,
              conv = ", bquote(.(res$convergence)), ": see optim() for details")
    
    it <- res$counts[1]   #number of iterations
    lL <- - res$value     #maximum value of the marginal loglikelihood
      if (frailty == "possta")
        lL <- lL + correct * log(10) * obsdata$ncl
  }
  
  #----- Recover the estimates ------------------------------------------------#
    #heterogeneity parameter
  if (frailty %in% c("gamma", "ingau")) {
    theta <- exp(res$par[1])
    nu <- NULL
  }
  else if (frailty == "possta") {
    nu <- exp(-exp(res$par[1]))
    theta <- NULL
  }
  else if (frailty == "none"){
    theta <- NULL
    nu <- NULL
  }
    #baseline hazard parameter(s)
  if (dist == "exponential") {
    lambda <- exp(res$par[nFpar+1])
    ESTIMATE <- c(lambda=lambda)
  }
  else if (dist == "weibull") {
    rho <- exp(res$par[nFpar+1])
    lambda <- exp(res$par[nFpar+2])
    ESTIMATE <- c(rho=rho, lambda=lambda)
  }
  else if (dist == "gompertz") {
    gamma <- exp(res$par[nFpar+1])
    lambda <- exp(res$par[nFpar+2])
    ESTIMATE <- c(gamma=gamma, lambda=lambda)
  }
  else if (dist == "lognormal") {
    mu <- res$par[nFpar+1]
    sigma <- exp(res$par[nFpar+2])
    ESTIMATE <- c(mu=mu, sigma=sigma)
  }
  else if (dist == "loglogistic") {
    alpha <- res$par[nFpar+1]
    kappa <- exp(res$par[nFpar+2])
    ESTIMATE <- c(alpha=alpha, kappa=kappa)
  }
    #regression parameter(s)
  if (nRpar == 0)
    beta <- NULL
  else{
    beta <- res$par[-(1:(nFpar + nBpar))]
    names(beta) <- paste("beta", names(obsdata$x), sep=".")[-1]
  }
    #all together
  ESTIMATE <- c(theta=theta,
                nu=nu,
                ESTIMATE,
                beta=beta)
  
  #----- Recover the standard errors ------------------------------------------#
  if ((frailty == "none") && is.null(inip)) {
    var <- coxMod$var
    
    if (nRpar == 0)
      seBeta <- NULL
    else {
      seBeta <- sqrt(var[1:nRpar])
      PVAL <- c(rep(NA, nFpar+nBpar), 
                2 * pnorm(q=- abs(beta / seBeta)))
    }
      
    if (dist == "exponential") {
      seLambda <- deltamethod(g=~exp(- x1),
                              mean=logscale,
                              cov=var["log(scale)", "log(scale)"],
                              ses=TRUE)
      STDERR <- c(seLambda=seLambda)
    }
    else if (dist == "weibull") {
      seRho <- deltamethod(g=~exp(x1), 
                           mean=logshape, 
                           cov=var["log(shape)", "log(shape)"], 
                           ses=TRUE)
      seLambda <- deltamethod(g=~exp(- exp(x2) * x1),
                              mean=c(logscale, logshape),
                              cov=var[c("log(scale)","log(shape)"),
                                      c("log(scale)","log(shape)")],
                              ses=TRUE)
      STDERR <- c(seRho=seRho, seLambda=seLambda)
    }
    else if (dist == "gompertz") {
      seGamma <- deltamethod(g=~exp(- x1),
                              mean=logscale,
                              cov=var["log(scale)", "log(scale)"],
                              ses=TRUE)
      seLambda <- deltamethod(g=~exp(x1),
                              mean=logshape,
                              cov=var["log(shape)", "log(shape)"],
                              ses=TRUE)
      STDERR <- c(seGamma=seGamma, seLambda=seLambda)
    }
    else if (dist == "lognormal") {
      seMu <- sqrt(var["log(scale)", "log(scale)"])
      seSigma <- deltamethod(g=~exp(- x1), 
                           mean=logshape, 
                           cov=var["log(shape)", "log(shape)"], 
                           ses=TRUE)
      STDERR <- c(seMu=seMu, seSigma=seSigma)
    }
    else if (dist == "loglogistic") {
      seAlpha <- deltamethod(g=~- exp(x2) * x1,
                              mean=c(logscale, logshape),
                              cov=var[c("log(scale)","log(shape)"),
                                      c("log(scale)","log(shape)")],
                              ses=TRUE)
      seKappa <- deltamethod(g=~exp(x1), 
                           mean=logshape, 
                           cov=var["log(shape)", "log(shape)"], 
                           ses=TRUE)
      STDERR <- c(seAlpha=seAlpha, seKappa=seKappa)
    }
    STDERR <- c(STDERR,
                se.beta=seBeta)
  }
  else {
    var <- try(diag(solve(res$hessian)), silent=TRUE)
    if (class(var) == "try-error") {
      warning(var[1])
      STDERR <- rep(NA, nFpar + nBpar + nRpar)
      PVAL <- rep(NA, nFpar + nBpar + nRpar)
    }
    else {
      if (any(var <= 0))
        warning(paste("negative variances have been replaced by NAs\n",
                      "Please, try other initial values",
                      "or another optimisation method"))
        #heterogeneity parameter(s)
      if (frailty %in% c("gamma", "ingau")) {
        seTheta <- ifelse(var[1] > 0,
                          sqrt(var[1] * theta^2),
                          NA)
        seNu <- NULL
      }
      else if (frailty == "possta") {
        seNu <- ifelse(var[1] > 0,
                       sqrt(var[1] * (nu * log(nu))^2),
                       NA)
        seTheta <- NULL
      }
        #baseline hazard parameter(s)
      if (dist == "exponential") {
        seLambda <- ifelse(var[nFpar+1] > 0,
                           sqrt(var[nFpar+1] * lambda^2),
                           NA)
        STDERR <- c(seLambda=seLambda)
      } 
      else if (dist == "weibull") {
        seRho <- ifelse(var[nFpar+1] > 0,
                        sqrt(var[nFpar+1] * rho^2),
                        NA)
        seLambda <- ifelse(var[nFpar+2] > 0,
                           sqrt(var[nFpar+2] * lambda^2),
                           NA)
        STDERR <- c(seRho=seRho, seLambda=seLambda)
      }
      else if (dist == "gompertz") {
        seGamma <- ifelse(var[nFpar+1] > 0,
                          sqrt(var[nFpar+1] * gamma^2),
                          NA)
        seLambda <- ifelse(var[nFpar+2] > 0,
                           sqrt(var[nFpar+2] * lambda^2),
                           NA)
        STDERR <- c(seGamma=seGamma, seLambda=seLambda)
      }
      else if (dist == "lognormal") {
        seMu <- ifelse(var[nFpar+1] > 0,
                       sqrt(var[nFpar+1]),
                       NA)
        seSigma <- ifelse(var[nFpar+2] > 0,
                          sqrt(var[nFpar+2] * sigma^2),
                          NA)
        STDERR <- c(seMu=seMu, seSigma=seSigma)
      }
      else if (dist == "loglogistic") {
        seAlpha <- ifelse(var[nFpar+1] > 0,
                          sqrt(var[nFpar+1]),
                          NA)
        seKappa <- ifelse(var[nFpar+2] > 0,
                          sqrt(var[nFpar+2] * kappa^2),
                          NA)
        STDERR <- c(seAlpha=seAlpha, seKappa=seKappa)
      }
        #regression parameter(s)
      if (nRpar == 0)
        seBeta <- NULL
      else {
        seBeta <- numeric(nRpar)
        varBeta <- var[-(1:(nFpar + nBpar))]
        for(i in 1:nRpar) {
          seBeta[i] <- ifelse(varBeta[i] > 0,
                              sqrt(varBeta[i]),
                              NA)
        }
        PVAL <- c(rep(NA, nFpar+nBpar), 
                  2 * pnorm(q=- abs(beta / seBeta)))
      }
    
      
        #all together
      STDERR <- c(STDERR,
                  se.beta=seBeta)
      if (frailty != "none")
        STDERR <- c(se.theta=seTheta,
                    se.nu=seNu,
                    STDERR)
        
    }
  }
        
  #----- Output ---------------------------------------------------------------#
  resmodel <- cbind(ESTIMATE=ESTIMATE, SE=STDERR)
  rownames(resmodel) <- gsub("beta.","", rownames(resmodel))
  
  if (nRpar > 0)
    resmodel <- cbind(resmodel, "p-val"= PVAL)

  class(resmodel) <- c("parfm", class(resmodel))
  attributes(resmodel) <- c(attributes(resmodel), list(
    convergence = res$convergence,
    it          = it,
    extime      = extime,
    nobs        = nrow(data),
    shared      = (nrow(data) > obsdata$ncl),
    loglik      = lL,
    dist        = dist,
    cumhaz      = attributes(Mloglikelihood(p=res$par,
                                 obs=obsdata, dist=dist, frailty=frailty,
                                 correct=correct))$cumhaz,
    di          = obsdata$di,
    frailty     = frailty,
    clustname   = cluster,
    correct     = correct))
  if (frailty != "none") {
    names(attr(resmodel, "cumhaz")) <-
      names(attr(resmodel, "di")) <- unique(obsdata$cluster)
  }
  
  if (showtime){
    cat("\nExecution time:", extime, "second(s) \n")
  }
  
  return(resmodel)
}

