\name{ecf}
\docType{data}
\alias{ecf}
\title{East Coast Fever transmission dynamics}
\description{
\emph{Theileriosis} or East Coast Fever (ECF) is a major cattle disease 
  in East and Southern Africa. 
The disease is caused by \emph{Theileria parva} which is transmitted by the ticks
\emph{Rhipicephalus appendiculatus} and the closely related
\emph{Rhipicephalus zambeziensis}.

In order to study the transmission of the disease, cows were followed up 
  from birth until time of first ECF contact in Nteme, 
  a region in Southern Zambia. 
On weekly basis, blood is collected and tested for the presence of antibodies 
  to \emph{Theileria parva} using the Indirect Fluorescent Antibody test. 
After three consecutive positive IFA test results an animal 
  is considered seroconverted.
The time to first ECF contact is defined as the timespan from birth 
  to one month before the first of the three consecutive positive test results 
  (i.e., we set the date of contact with \emph{Theileria parva} one month before 
  the time of the first positive test result). 
Animals are followed up until one year after birth; if they do not seroconvert 
  by that time, their time of seroconversion is right-censored. 
We also consider the binary covariate breed. 
}
\usage{data(ecf)}
 \format{A dataframe containing 212 observations.
 \describe{
    \item{Cowid:}{Cow's identifyier.}
    \item{Time:}{Time to ECF contact (in days).}
    \item{Status:}{Censored (0) or observed (1) event time.}
    \item{Breed:}{The cow's breed.}
  }
}
\source{Example 1.1 of Duchateau an Janssen (2008)
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
\donttest{
library(parfm)
data(ecf)
head(ecf)

################################################################################
#Example 3.8: Population and conditional hazard for time to ECF contact        #
#Duchateau and Janssen (2008, page 113)                                        #
################################################################################
pfm1 <- parfm(Surv(Time, Status) ~ 1, cluster = "Cowid", data = ecf,
              dist = "weibull", frailty = "gamma")
pfm2 <- parfm(Surv(Time, Status) ~ Breed, cluster = "Cowid", data = ecf,
              dist = "weibull", frailty = "gamma")
curve(pfm1["lambda", 1] * pfm1["rho", 1] * x ^ (pfm1["rho", 1] - 1),
      from = 0, to = 400, ylim = c(0, .15), 
      ylab = "Hazard function", xlab = "Time (days)")
curve(qgamma(.75, shape = 1 / pfm1["theta", 1],
             scale = pfm1["theta", 1]) * pfm1["lambda", 1] * pfm1["rho", 1] *
                x ^ (pfm1["rho", 1] - 1),
      add = TRUE, lty = 2)
curve(qgamma(.25, shape=1 / pfm1["theta", 1],
             scale = pfm1["theta", 1]) * pfm1["lambda", 1] * pfm1["rho", 1] *
               x ^ (pfm1["rho", 1] - 1),
      add = TRUE, lty = 3)
curve(pfm1["lambda", 1] * pfm1["rho", 1] *
        x ^ (pfm1["rho", 1] - 1) / (
            1 + pfm1["theta", 1] * pfm1["lambda", 1] * x ^ (pfm1["rho", 1])
      ),
      add = TRUE, lwd = 2)
legend("top", lwd = c(1, 1, 1, 2), lty = c(1, 2, 3, 1), ncol = 2,
       legend = c("Mean frailty", "Q75 frailty", "Q25 frailty", "Population"))
}}
\keyword{datasets}
\keyword{ecf}
\keyword{cows}
\keyword{Theileriosis}
\keyword{East Coast Fever}
\keyword{infection}
\keyword{frailty}