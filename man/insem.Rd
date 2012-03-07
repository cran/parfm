\name{insem}
\docType{data}
\alias{insem}
\title{Time to first insemination in dairy heifer cows 
  without time varying covariates}
\description{
In a dairy farm, the calving interval (the time between two calvings) 
  should be optimally between 12 and 13 months. 
One of the main factors determining the length of the calving interval 
  is the time from parturition to the time of first insemination. 

The objective of this study is to look for cow factors that might predict 
  the time to first insemination, 
  so that actions can be taken based on these predictors. 
As no inseminations take place in the first 29 days after calving, 
  we subtract 29 days (and not 30 days as the first event would then have 
  first insemination time zero) from the time to first insemination 
  since at risk time starts only then. 
Cows which are culled without being inseminated are censored 
  at their culling time. 
Furthermore, cows that are not yet inseminated 300 days after calving 
  are censored at that time. 

The time to first insemination is studied in dairy cows 
  as a function of covariates that are fixed over time.   
An example ofsuch a covariate is the parity of the cow, 
  corresponding to the number of calvings the cow has had already. 
As we observe only one lactation period for each cow in the study, 
  it is indeed a constant cow characteristic 
  within the time framework of the study. 

We dichotomise parity into primiparous cows or heifers
  (only one calving (\code{heifer=1})) and multiparous cows
  (more than one calving (\code{heifer=0})). 
Other covariates that are used in the analysis are 
  the different milk constituents such as protein and 
  ureum concentration at parturition. 
}
\usage{data(insem)}
 \format{A dataframe containing 10513 observations.
 \describe{
    \item{Cowid:}{Cow's identifyier.}
    \item{Time:}{Time to first insemination (in days).}
    \item{Status:}{Censored (0) or observed (1) event time.}
    \item{Herd:}{The herd to which the cow belongs.}
    \item{Urem:}{Milk urem concentration (\%)
      at the start of the lactation period.}
    \item{Protein:}{Protein concentration (\%)
      at the start of the lactation period.}
    \item{Parity:}{The number of calvings.}
    \item{Heifer:}{Multiparous cow (0) or primiparous cow (1).}
  }
}
\source{Example 1.8 of Duchateau an Janssen (2008)
  \url{http://www.vetstat.ugent.be/research/frailty/datasets/}}
\note{
  These data are downloaded from
    \url{http://www.vetstat.ugent.be/research/frailty/datasets/}.
  They are simulated data, with exactly
    the same structure as the real data used in the book,
    that could not be made publicly available.    
}
\references{
  Duchateau L., Janssen P. (2008). \emph{The frailty model. Springer}.
  New York: Springer--Verlag.
}
\examples{
\donttest{
library(parfm)
data(insem)
head(insem)

insem$TimeMonths <- insem$Time * 12 / 365.25
  
################################################################################
#Example 2.1: The parametric proportional hazards frailty model for the time   #
#to first insemination based on marginal likelihood maximisation               #
#Duchateau and Janssen (2008, page 50)                                         #
################################################################################
pfm <- parfm(Surv(TimeMonths, Status)~Heifer, cluster="Herd", data=insem,
             dist="weibull", frailty="gamma")
pfm

par(mfrow=c(2, 2))

### - Hazard functions - ###
# multiparous cows
curve((365.25 / 12)^(-pfm["rho", 1]) *
      pfm["lambda", 1] * pfm["rho", 1] * x^(pfm["rho", 1] - 1),
      from=0, to=400, ylim=c(0, .14), 
      main="Multiparous cows",
      ylab="Hazard function", xlab="Time to first insemination (days)")
curve(qgamma(.95, shape=1/pfm["theta", 1], scale=pfm["theta", 1]) *
      (365.25 / 12)^(-pfm["rho", 1]) *
      pfm["lambda", 1] * pfm["rho", 1] * x^(pfm["rho", 1] - 1),
      add=TRUE, lty=4)
curve(qgamma(.05, shape=1/pfm["theta", 1], scale=pfm["theta", 1]) *
      (365.25 / 12)^(-pfm["rho", 1]) *
      pfm["lambda", 1] * pfm["rho", 1] * x^(pfm["rho", 1] - 1),
      add=TRUE, lty=4)
      
# primiparous cows
curve(exp(pfm["Heifer", 1]) * 
      (365.25 / 12)^(-pfm["rho", 1]) *
      pfm["lambda", 1] * pfm["rho", 1] * x^(pfm["rho", 1] - 1),
      from=0, to=400, ylim=c(0, .14), 
      main="Primiparous cows", 
      ylab="Hazard function", xlab="Time to first insemination (days)")
curve(qgamma(.95, shape=1/pfm["theta", 1], scale=pfm["theta", 1]) *
      exp(pfm["Heifer", 1]) * 
      (365.25 / 12)^(-pfm["rho", 1]) *
      pfm["lambda", 1] * pfm["rho", 1] * x^(pfm["rho", 1] - 1),
      add=TRUE, lty=4)
curve(qgamma(.05, shape=1/pfm["theta", 1], scale=pfm["theta", 1]) *
      exp(pfm["Heifer", 1]) * 
      (365.25 / 12)^(-pfm["rho", 1]) *
      pfm["lambda", 1] * pfm["rho", 1] * x^(pfm["rho", 1] - 1),
      add=TRUE, lty=4)
       
    
### - Cumulative distribution functions - ###
# multiparous cows
curve(1 - exp(-
        (365.25 / 12)^(-pfm["rho", 1]) *
        pfm["lambda", 1] * x^(pfm["rho", 1])),
      from=0, to=400, ylim=c(0, 1), 
      main="Multiparous cows", 
      ylab="Cumulative distribution function", 
      xlab="Time to first insemination (days)")
curve(1 - exp(-
        qgamma(.95, shape=1/pfm["theta", 1], scale=pfm["theta", 1]) *
        (365.25 / 12)^(-pfm["rho", 1]) *
        pfm["lambda", 1] * x^(pfm["rho", 1])),
      add=TRUE, lty=4)
curve(1 - exp(-
        qgamma(.05, shape=1/pfm["theta", 1], scale=pfm["theta", 1]) *
        (365.25 / 12)^(-pfm["rho", 1]) *
        pfm["lambda", 1] * x^(pfm["rho", 1])),
      add=TRUE, lty=4)
      
# primiparous cows
curve(1 - exp(-
        exp(pfm["Heifer", 1]) * 
        (365.25 / 12)^(-pfm["rho", 1]) *
        pfm["lambda", 1] * x^(pfm["rho", 1])),
      from=0, to=400, ylim=c(0, 1), 
      main="Primiparous cows", 
      ylab="Cumulative distribution function", 
      xlab="Time to first insemination (days)")
curve(1 - exp(-
        qgamma(.95, shape=1/pfm["theta", 1], scale=pfm["theta", 1]) *
        exp(pfm["Heifer", 1]) * 
        (365.25 / 12)^(-pfm["rho", 1]) *
        pfm["lambda", 1] * x^(pfm["rho", 1])),
      add=TRUE, lty=4)
curve(1 - exp(-
        qgamma(.05, shape=1/pfm["theta", 1], scale=pfm["theta", 1]) *
        exp(pfm["Heifer", 1]) * 
        (365.25 / 12)^(-pfm["rho", 1]) *
        pfm["lambda", 1] * x^(pfm["rho", 1])),
      add=TRUE, lty=4)
      
      
### - Density of the median time - ###
fM <- function(x, heifer=0) {
  RHO <- pfm["rho", 1]
  LAMBDAd <- (365.25 / 12)^(-pfm["rho", 1]) * pfm["lambda", 1]
  THETA <- pfm["theta", 1]
  if (heifer) {
    eBETA <- exp(pfm["Heifer", 1])
  } else eBETA <- 1  
  
  RHO * (log(2) / (THETA * LAMBDAd * eBETA))^(1 / THETA) *
  x^(-1 - RHO / THETA) *
  exp(-log(2) / (THETA * x^RHO * LAMBDAd * eBETA)) /
  gamma(1 / THETA)
}

par(mfrow=c(1, 1))
curve(fM, 0, 300, xlab="Median time to first insemination (days)",
      ylab="Density function of the median")
curve(fM(x, heifer=1), add=TRUE, lty=3)
legend("topright", legend=c("Multiparous", "Primiparous"),
       lty=c(1, 3), bty="n")
}
}
\keyword{datasets}
\keyword{insem}
\keyword{cows}
\keyword{insemination}
\keyword{frailty}