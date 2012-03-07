\name{insemtvc}
\docType{data}
\alias{insemtvc}
\title{Time to first insemtvcination in dairy heifer cows 
  with time varying covariates}
\description{
In a dairy farm, the calving interval (the time between two calvings) 
  should be optimally between 12 and 13 months. 
One of the main factors determining the length of the calving interval 
  is the time from parturition to the time of first insemtvcination. 

The objective of this study is to look for cow factors that might predict 
  the time to first insemtvcination, 
  so that actions can be taken based on these predictors. 
As no insemtvcinations take place in the first 29 days after calving, 
  we subtract 29 days (and not 30 days as the first event would then have 
  first insemtvcination time zero) from the time to first insemtvcination 
  since at risk time starts only then. 
Cows which are culled without being inseminated are censored 
  at their culling time. 
Furthermore, cows that are not yet inseminated 300 days after calving 
  are censored at that time. 

The time to first insemtvcination is studied in dairy cows 
  as a function of time-varying covariates.
Examples of such covariates are the protein and ureum concentration. 
The protein and ureum concentrations are measured, during the experiment, 
  at a number of points in time. 
It might be more adequate to model the hazard at a particular time 
  using the concentration at that particular point in time. 
In order to accommodate for time-varying covariates, 
  the risk time for each cow is split into time intervals 
  with a start and an end time and in each such interval 
  a constant value for the concentration. 
We therefore have for each cow as many data lines as there are risk intervals. 
}
\usage{data(insemtvc)}
 \format{A dataframe containing 10513 observations.
 \describe{
    \item{Cowid:}{Cow's identifyier.}
    \item{Start:}{Start time of the interval (in days).}
    \item{End:}{End time of the interval (in days).}
    \item{Status:}{Censored (0) or observed (1) event time.}
    \item{Herd:}{The herd to which the cow belongs.}
    \item{Urem:}{Milk urem concentration (\%)
      at the begin time of the interval.}
    \item{Protein:}{Protein concentration (\%)
      at the begin time of the interval.}
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
\keyword{datasets}
\keyword{insemtvc}
\keyword{cows}
\keyword{insemtvcination}
\keyword{frailty}
