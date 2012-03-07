\name{culling}
\docType{data}
\alias{culling}
\title{Culling of dairy heifer cows}
\description{
  The time to culling is studied in heifers as a function 
  of the somatic cell count (SCC) measured between 5 and 15 days
 (measurement day) after calving. 
 
High somatic cell count (we use the logarithm of somatic cell count as covariate)
  might be a surrogate marker for intramammary infections. 
Heifers which have intramammary infections or which are expected to
  develop intramammary infections in the future are quite expensive to keep
  due to the high costs for drugs and the loss in milk production.
  
Cows are followed up for an entire lactation period (roughly 300--50 days) 
  and, if they are still alive at the end of the lactation period,
  they are censored at that time. 
Cows are further clustered within herds and this clustering needs
 to be taken into account as culling policy and also SCC
 in early lactation might differ substantially between the herds.  
}
\usage{data(culling)}
 \format{A dataframe containing 13836 observations.
 \describe{
    \item{Cowid:}{Cow's identifyier.}
    \item{Time:}{Time to culling (in days).}
    \item{Status:}{Censored (0) or observed (1) event time.}
    \item{Herd:}{Herd identifyier.}
    \item{Timeasses:}{SCC measurement day.}
    \item{LogSCC:}{Logarithm of the somatic cell count.}
  }
}
\source{Example 1.7 of Duchateau an Janssen (2008)
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
data(culling)
head(culling)
culling <- culling[culling$Time > 0,]
culling$TimeMonths <- culling$Time * 12 / 365.25

coxmod <- parfm(Surv(TimeMonths, Status)~LogSCC, data=culling,
                dist="exponential", frailty="none")
coxmod

pfmod <- parfm(Surv(TimeMonths, Status)~LogSCC, cluster="Herd", data=culling,
               dist="exponential", frailty="gamma")
pfmod
}
}
\keyword{datasets}
\keyword{culling}
\keyword{cows}
\keyword{frailty}

