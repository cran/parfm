\name{mastitis}
\docType{data}
\alias{mastitis}
\title{Correlated infection times in four cow udder quarters}
\description{
Mastitis, the infection of the udder, is economically 
  the most important disease in the dairy sector of the western world. 
Mastitis can be caused by many organisms, most of them bacteria, 
  such as \emph{Escherichia coli}, \emph{Streptococcus uberis}, 
  and \emph{Staphylococcus aureus}. 
Since each udder quarter is separated from the three other quarters, 
  one quarter might be infected with the other quarters free of infection. 
  
In an extensive study, 100 cows are followed up for infections. 
The objective of this observational study is to estimate the incidence 
  of the different organisms causing mastitis in the dairy cattle population 
  in Flanders. 
Also the correlation between the infection times of the four udder quarters of 
  a cow is an important parameter to take preventive measures against mastitis.
With high correlation, a lot of attention should be given to 
  the uninfected udder quarters of a cow that has an infected quarter. 
From each quarter, a milk sample is taken monthly and is screened 
  for the presence of different bacteria. 
  
We model the time to infection with any bacteria, with the cow being
  the cluster and the quarter the experimental unit within the cluster. 
Observations can be right-censored if no infection occurs before the end 
  of the lactation period, which is roughly 300-350 days but different 
  for every cow, or if the cow is lost to follow-up during the study, 
  for example due to culling. 
Due to the periodic follow-up, udder quarters that experience an event are
  interval-censored with lowerbound the time of the last milk sample
  with a negative result and upperbound the time of the first milk sample
  with a positive result.   
In some examples, the midpoint (average of lowerbound and upperbound 
  of the interval) is used for simplicity; 
  in other examples the interval-censored nature of the data
  is taken into account. 
  
In the analysis, two types of covariates are considered. 
Cow level covariates take the same value for every udder quarter of the cow
  (e.g., number of calvings or parity). 
Several studies have shown that prevalence as well as incidence of 
  intramammary infections increase with parity . 
Several hypotheses have been suggested to explain these findings, e.g., 
  teat end condition deteriorates with increasing parity. 
Because the teat end is a physical barrier that prevents organisms from 
  invading the udder, impaired teat ends make the udder more vulnerable 
  for intramammary infections. 
For simplicity, parity is dichotomised into primiparous cows (\code{heifer=1}) 
  and multiparous cows (\code{heifer=0}). 
Udder quarter level covariates change within the cow (e.g., position of the
  udder quarter, front or rear). 
The difference in teat end condition between front and rear quarters 
  has also been put forward to explain the difference in infection status.
}
\usage{data(mastitis)}
 \format{A dataframe containing 400 observations.
 \describe{
    \item{Cowid:}{Cow's identifyier.}
    \item{Lower:}{Lower bound of time to infection (in days).}
    \item{Upper:}{Upper bound of time to infection (in days).}
    \item{Midpoint:}{Midpoint between \code{Lower} and \code{Upper} (in days).}
    \item{Status:}{Censored (0) or observed (1) event time.}
    \item{Heifer:}{Multiparous cow (0) or primiparous cow (1).}
    \item{Quarter:}{Tthe udder quarter
      (\code{LF}=Left-Front,
      \code{LR}=Left-Rear,
      \code{RF}=Right-Front,
      \code{RR}=Right-Rear).}
  }
}
\source{Example 1.4 of Duchateau an Janssen (2008)
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
data(mastitis)
head(mastitis)

mastitis$timeto <- as.numeric((mastitis$Midpoint * 4 / 365.25))


################################################################################
# Example 4.4: The gamma frailty model for the udder quarter infection data    #
# Duchateau and Janssen (2008, page 136)                                       #
################################################################################
modParfm <- parfm(Surv(timeto, Status) ~ Heifer,
                  cluster = "Cowid",
                  dist = "weibull",
                  frailty = "gamma",
                  data = mastitis)
modParfm


################################################################################
# Example 4.7 The inverse Gaussian frailty model for the udder quarter         #
# infection data                                                               #
# Duchateau and Janssen (2008, page 156)                                       #
################################################################################
mastitis <- data.frame(mastitis, 
                       timeto = as.numeric((mastitis$Midpoint * 4 / 365.25)))
modParfm <- parfm(Surv(timeto, Status) ~ Heifer,
                  cluster = "Cowid",
                  dist = "weibull",
                  frailty = "ingau",
                  data = mastitis)
modParfm


################################################################################
# Example 4.10 The positive stable frailty model for the udder quarter         #
# infection data                                                               #
# Duchateau and Janssen (2008, page 169)                                       #
################################################################################
mastitis <- data.frame(mastitis, 
                       timeto = as.numeric((mastitis$Midpoint * 4 / 365.25)))
modParfm <- parfm(Surv(timeto, Status) ~ Heifer,
                  cluster = "Cowid",
                  dist = "weibull",
                  frailty = "possta",
                  data = mastitis)
modParfm
}}
\keyword{datasets}
\keyword{mastitis}
\keyword{cows}
\keyword{infection}
\keyword{udder}
\keyword{frailty}
