\name{reconstitution}
\docType{data}
\alias{reconstitution}
\title{Reconstitution of blood--milk barrier after reconstitution}
\description{
When an udder quarter of a cow is infected (reconstitution), 
  the blood-milk barrier is partially destroyed and particular ions 
  can flow freely from blood to milk and vice versa, leading to higher 
  concentrations of, for instance, the sodium concentration Na. 

The objective of this study is to demonstrate that the local application 
  of a drug based on corticosteroids decreases the time to reconstitution 
  of the blood-milk barrier in dairy cows. 
We therefore consider as outcome the time until the Na concentration goes
  below a certain threshold (a concentration below the threshold value
  is considered to be normal again). 
Each udder quarter is separated from the three other quarters so that a quarter
  can be used as experimental unit to which a treatment is assigned. 
The Na concentration in each of the experimental units is followed up. 
The rear udder quarters of 100 cows are experimentally infected 
  with \code{Escherichia coli}. 
After nine hours, one of the two infected udder quarters is treated locally 
  with the active compound whereas the other is treated with placebo. 
Cows are followed up for 6.5 days, and are censored at that point in time 
  if the Na concentration is still above the threshold level. 
  
We further include parity in the study as covariate. 
The parity of a cow is the number of calvings (and therefore the number
  of lactation periods) that the cow has already experienced. 
Parity is often converted into a binary covariate, grouping all the cows 
  with more than one calving in the group of multiparous cows (\code{heifer=0})
  compared to the group of primiparous cows or heifers, cows with 
  only one calving (\code{heifer=1}). 
}
\usage{data(reconstitution)}
 \format{A dataframe containing 200 observations.
 \describe{
    \item{Cowid:}{Cow's identifyier.}
    \item{Time:}{The time to reconstitution (in days).}
    \item{Status:}{Censored (0) or observed (1) event time.}
    \item{Drug:}{Drug or placebo application.}
    \item{Heifer:}{Multiparous cow (0) or primiparous cow (1).}
  }
}
\source{Example 1.3 of Duchateau an Janssen (2008)
  \url{http://www.vetstat.ugent.be/research/frailty/datasets/}}
\note{
  These data are downloaded from
    \url{http://www.vetstat.ugent.be/research/frailty/datasets/}.
  They are simulated data, with exactly
    the same structure as the real data used in the book,
    that could not be made publicly available.    
}
\references{
  Duchateau L, Janssen P (2008). \emph{The frailty model. Springer}.
  New York: Springer--Verlag.
}
\examples{
data(reconstitution)
head(reconstitution)
                

################################################################################
# Example 3.1: The [...] frailty model for the time to blood-milk barrier      #
# reconstitution: the effect of treatment                                      #
# Duchateau and Janssen (2008, page 79)                                        #
################################################################################
pfmDrug <- parfm(Surv(Time, Status) ~ Drug,
                 cluster = "Cowid", dist = "exponential",
                 frailty = "gamma", data = reconstitution)
ci.parfm(pfmDrug)
                

################################################################################
# Example 3.2: The [...] frailty model for the time to blood-milk barrier      #
# reconstitution: the heifer effect                                            #
# Duchateau and Janssen (2008, page 82)                                        #
################################################################################
pfmHeifer <- parfm(Surv(Time, Status) ~ Heifer,
                   cluster = "Cowid", dist = "exponential",
                   frailty = "gamma", data = reconstitution)
ci.parfm(pfmHeifer)
}
\keyword{datasets}


\keyword{reconstitution}
\keyword{reconstitution}
\keyword{cows}
\keyword{infection}
\keyword{udder}
\keyword{frailty}