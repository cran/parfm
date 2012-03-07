\name{asthma}
\docType{data}
\alias{asthma}
\title{Recurrent asthma attacks in children}
\description{
  Asthma is occurring more and more frequently in very young children
  (between 6 and 24 months).
  Therefore, a new application of an existing anti-allergic drug is administered
  to children who are at higher risk to develop asthma in order to prevent it.
  A prevention trial is set up with such children randomised to placebo or drug, 
  and the asthma events that developed over time are recorded in a diary. 
  Typically, a patient has more than one asthma event. 
  The different events are thus clustered within a patient and are ordered 
  in time. 
  This ordering can be taken into account in the model.
  Such data can be presented in different formats, but here, we choose to use
  the calendar time representation. In the calendar time representation, 
  the time at risk for a particular event is the time from the end
  of the previous event (asthma attack) to the start of the next event
  (start of the next asthma attack). 
  In describing recurrent event data, we need a somewhat more complex 
  data structure to keep track of the sequence of events within a patient.
  A particular patient has different periods at risk during the total 
  observation period which are separated either by an asthmatic event that 
  lasts one or more days or by a period in which the patient was not 
  under observation. 
  The start and end of each such risk period is required, together with the 
  status indicator to denote whether the end of the risk period corresponds to
  an asthma attack or not.
}
\usage{data(asthma)}
 \format{A dataframe containing 1776 observations.
 \describe{
    \item{Patid:}{Patient's identifyier.}
    \item{Begin:}{Time of end of the previous asthma attack (in days).}
    \item{End:}{Asthma attack or censoring time. (in days)}
    \item{Status:}{Censored (0) or observed (1) event time.}
    \item{Drug:}{placebo (0) or drug (1).}
    \item{Fevent:}{First observation of the patient? 1=yes, 0=no.}
  }
}
\author{
  F. Rotolo and M. Munda. 
  Original text and data by L. Duchateau and P. Janssen.
}
\source{Example 1.9 of Duchateau an Janssen (2008)
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
data(asthma)
head(asthma)
asthma <- asthma[asthma$Fevent==0,]

################################################################################
#Example 2.4: The frailty model with the Weibull baseline for the recurrent    #
#asthma data based on marginal likelihood maximisation                         #
#Duchateau and Janssen (2008, page 56)                                         #
################################################################################
# Calendar time
parfm(Surv(Begin, End, Status)~Drug, cluster="Patid", 
      data=asthma,
      dist="weibull", frailty="gamma")

# Gap time
asthma$time <- asthma$End - asthma$Begin
parfm(Surv(time, Status)~Drug, cluster="Patid", data=asthma,
      dist="weibull", frailty="gamma")
}
}
      
\keyword{datasets}
\keyword{asthma}
\keyword{frailty}

