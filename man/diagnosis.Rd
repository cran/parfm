\name{diagnosis}
\docType{data}
\alias{diagnosis}
\title{Diagnosis of fracture healing}
\description{
Medical imaging has become an important tool in the veterinary hospital
  to assess whether and when a fracture has healed. 
The standard technique in dogs is based on radiography (RX). 
Newer techniques based on ultrasound (US) are cheaper 
  and do not require radioprotection. 

To investigate the performance of US for this purpose and to compare it to RX,
  Risselada et al. (2006) set up a trial in which fracture healing is evaluated 
  by both US and RX. 
In total, 106 dogs, treated in the veterinary university hospital of Ghent, 
  are included in the trial and evaluated for time to fracture healing 
  with the two techniques. 
Only 7 dogs are censored for time to fracture healing evaluated by RX; 
  no censoring occurs for time to fracture healing evaluated by US. 
The censoring is due to the fact that dog owners do not show up anymore. 
}
\usage{data(diagnosis)}
 \format{A dataframe containing 212 observations.
 \describe{
    \item{Dogid:}{Dog's identifyier.}
    \item{Time:}{Time to diagnosis (in days).}
    \item{Status:}{Censored (0) or observed (1) event time.}
    \item{Method:}{Diagnostic technique: either
      \code{RX}, radiography, or \code{US}, ultrasound.}
  }
}
\source{Example 1.2 of Duchateau an Janssen (2008)
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
  
  Risselada, M., van Bree, H., Kramer, M., Chiers, K., Duchateau, L. 
  and Verleyen, P. (2006). 
  Evaluation of nonunion fractures in dogs by use of Bmode ultrasonography, 
  power Doppler ultrasonography, radiography, and histologic examination. 
  \emph{Am. J. Vet. Res.} 67, 1354--1361.
}
\examples{
library(parfm)
data(diagnosis)
head(diagnosis)

diagnosis$TimeMonths <- diagnosis$Time * 12 / 365.25

################################################################################
#Example 3.6: Shared gamma frailty models [...] for time to diagnosis          #
#of being healed                                                               #
#Duchateau and Janssen (2008, page 101)                                        #
################################################################################
WeiGam <- parfm(Surv(TimeMonths, Status)~Method, cluster="Dogid", data=diagnosis,
               dist="weibull", frailty="gamma")
\donttest{
WeiGam

curve(WeiGam["lambda", 1] * WeiGam["rho", 1] * x^(WeiGam["rho", 1] - 1),
      from=0, to=4, ylim=c(0, 25), 
      ylab="Hazard function", xlab="Time (months)")
curve(exp(WeiGam["Method", 1]) *
      WeiGam["lambda", 1] * WeiGam["rho", 1] * x^(WeiGam["rho", 1] - 1),
      add=TRUE, lty=2)
legend("topleft", lty=1:2, legend=c("US", "RX"))
      
################################################################################
#Example 4.8: Inverse Gaussian frailty models [...] for time to diagnosis      #
#of being healed                                                               #
#Duchateau and Janssen (2008, page 160)                                        #
################################################################################
}
WeiIG <- parfm(Surv(TimeMonths, Status)~Method, cluster="Dogid", data=diagnosis,
               dist="weibull", frailty="ingau")
\donttest{
WeiIG

curve(WeiIG["lambda", 1] * WeiIG["rho", 1] * x^(WeiIG["rho", 1] - 1),
      from=0, to=4, ylim=c(0, 25), 
      ylab="Hazard function", xlab="Time (months)")
curve(exp(WeiIG["Method", 1]) *
      WeiIG["lambda", 1] * WeiIG["rho", 1] * x^(WeiIG["rho", 1] - 1),
      add=TRUE, lty=2)
legend("topleft", lty=1:2, legend=c("US", "RX"))


################################################################################
#Example 4.11: Positive Stable frailty models [...] for time to diagnosis      #
#of being healed                                                               #
#Duchateau and Janssen (2008, page 172)                                        #
################################################################################

WeiPS <- parfm(Surv(TimeMonths, Status)~Method, cluster="Dogid", data=diagnosis,
               dist="weibull", frailty="possta")
WeiPS

curve(WeiPS["lambda", 1] * WeiPS["rho", 1] * x^(WeiPS["rho", 1] - 1),
      from=0, to=4, ylim=c(0, 25), 
      ylab="Hazard function", xlab="Time (months)")
curve(exp(WeiPS["Method", 1]) *
      WeiPS["lambda", 1] * WeiPS["rho", 1] * x^(WeiPS["rho", 1] - 1),
      add=TRUE, lty=2)
legend("topleft", lty=1:2, legend=c("US", "RX"))
}
}
\keyword{datasets}
\keyword{diagnosis}
\keyword{fracture}
\keyword{dogs}
\keyword{RX}
\keyword{radiography}
\keyword{US}
\keyword{ultrasound}
\keyword{frailty}