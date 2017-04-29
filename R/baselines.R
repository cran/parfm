################################################################################
#  Baseline hazard distributions                                               #
################################################################################
#                                                                              #
#  These are functions with parameters                                         #
#   - pars: the vector of parameters                                           #
#   - t   : the time point                                                     #
#   - what: the quantity to be returned by the function                        #
#           either  "H", the cumulated hazard,                                 #
#               or "lh", the log-hazard                                        #
#                                                                              #
#                                                                              #
#   Date:                 December 19, 2011                                    #
#   Last modification on: January 25, 2017                                     #
################################################################################


################################################################################
################################################################################


################################################################################
#                                                                              #
#   Weibull baseline hazard function                                           #
#                                                                              #
#   Parameters:                                                                #
#    [1] rho    > 0                                                            #
#    [2] lambda > 0                                                            #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \rho \lambda t ^ (\rho-1)                                          #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
################################################################################

weibull <- function(pars,
                    t, 
                    what) {
    if (what == "H")
        return(pars[2] * t ^ (pars[1]))
    else if (what == "lh")
        return(log(pars[1]) + log(pars[2]) + ((pars[1] - 1) * log(t)))
}


################################################################################
#                                                                              #
#   Fréchet or Inverse Weibull baseline hazard function                        #
#                                                                              #
#   Parameters:                                                                #
#    [1] rho    > 0                                                            #
#    [2] lambda > 0                                                            #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \lambda \rho t ^ -(\rho + 1) / (exp(\lambda t ^ -\rho) - 1)        #
#    H(t) = -log[1 - exp{- \lambda t ^ -\rho}]                                 #
#                                                                              #
#                                                                              #
#   Date:              June, 26, 2012                                          #
#   Last modification: January 31, 2017                                        #
################################################################################
inweibull <- frechet <- function(pars,
                                 t, 
                                 what) {
    if (what == "H")
        return(-log(1 - exp(-pars[2] * (t ^ -pars[1]))))
    else if (what == "lh")
        return(sum(log(pars[1:2])) - log(t) * (pars[1] + 1) -
                   log(exp(pars[2] * (t ^ -pars[1])) - 1))
}


################################################################################
#                                                                              #
#   Exponential baseline hazard function                                       #
#                                                                              #
#   Parameters:                                                                #
#    [1] lambda > 0                                                            #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \lambda                                                            #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
################################################################################

exponential <- function(pars,
                        t, 
                        what) {
    if (what == "H")
        return(pars * t)
    else if (what == "lh") 
        return(log(pars))
}



################################################################################
#                                                                              #
#   Gompertz baseline hazard function                                          #
#                                                                              #
#   Parameters:                                                                #
#    [1] gamma  > 0                                                            #
#    [2] lambda > 0                                                            #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \lambda exp(\gamma t)                                              #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
################################################################################

gompertz <- function(pars,
                     t, 
                     what) {
    if (what == "H") 
        return(pars[2] / pars[1] * (exp(pars[1] * t) - 1))
    else if (what == "lh") 
        return(log(pars[2]) + pars[1] * t)
}



################################################################################
#                                                                              #
#   Lognormal baseline hazard function                                         #
#                                                                              #
#   Parameters:                                                                #
#    [1] mu     \in \mathbb R                                                  #
#    [2] sigma  > 0                                                            #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \frac{ \phi(\log t -\mu \over \sigma) }{                           #
#             \sigma t [1-\Phi(\log t -\mu \over \sigma)] }                    #
#                                                                              #
#   with \phi the density and \Phi the distribution function                   #
#   of a standard Normal                                                       #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
################################################################################

lognormal <- function(pars,
                      t, 
                      what) {
    if (what == "H")  #return - log (S)
        return(- log(1 - plnorm(t, meanlog = pars[1], sdlog = pars[2])))
    else if (what == "lh")  #return log(f) - log(S)
        return(dlnorm(t, meanlog = pars[1], sdlog = pars[2], log = TRUE) -
                   log(1 - plnorm(t, meanlog = pars[1], sdlog = pars[2])))
}



################################################################################
#                                                                              #
#   Loglogistic baseline hazard function                                       #
#                                                                              #
#   Parameters:                                                                #
#    [1] alpha \in \mathbb R                                                   #
#    [2] kappa > 0                                                             #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \frac{ \exp(\alpha) \kappa t^{\kappa-1} }{                         #
#             1 + \exp(\alpha) t^\kappa }                                      #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
################################################################################

loglogistic <- function(pars,
                        t, 
                        what) {
    if (what == "H") 
        return(log(1 + exp(pars[1]) * t ^ (pars[2])))
    else if (what == "lh") 
        return(pars[1] + log(pars[2]) + (pars[2] - 1) * log(t) - 
                   log(1 + exp(pars[1]) * t ^ pars[2]) )
}




################################################################################
#                                                                              #
#   Log-skewNormal baseline hazard function                                    #
#   Azzalini, Adelchi. (1985). A class of distributions which includes         #
#     the normal ones. Scandinavian Journal of Statistics, 12:171-178          # 
#                                                                              #
#   Parameters:                                                                #
#    [1] xi \in \mathbb R, the location parameter                              #
#    [2] omega > 0, the scale parameter (called omega in original paper)       #
#    [3] alpha \in \mathbb R, the shape parameter                              #
#                                                                              #
#   Density:                                                                   #
#    f(t) = 2  / (\omega t) \phi_d((\log(t) - \xi) / \omega)                   #
#                                        \Phi(\alpha (log(t) - \xi) / \omega)  #
#                                                                              #
#                                                                              # 
#   Author: Andrea Callegaro                                                   #
#   Date: November, 22, 2016                                                   #
################################################################################

logskewnormal <- function(pars,
                          t, 
                          what) {
    # library(sn)
    #log-skew normal density
    dlsn <- Vectorize(function(t, pars) {
        dsn(log(t), xi = pars[1], omega = pars[2], alpha = pars[3]
        ) / t
    }, 't')
    #log-skew normal cdf
    plsn <- Vectorize(function(t, pars) {
        psn(log(t), xi = pars[1], omega = pars[2], alpha = pars[3])
    }, 't')
    
    if (what == "H")  #return -log(S)
        return(-log(1 - plsn(t, pars = pars)))
    
    else if (what == "lh")  #return (log(f) - log(S))
        return(log(dlsn(t, pars = pars)) -
                   log(1 - plsn(t, pars = pars)))
}


# setGeneric('basehaz', package = 'survival')
# setClass('parfm')
# setMethod('basehaz', signature(fit = 'parfm'), function(fit, centered = TRUE) {
#     return(fit)
# })
