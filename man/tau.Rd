\name{tau}
\alias{tau}
\title{Kendall's Tau for Parametric Frailty Models}
\description{Computes Kendall's Tau for Parametric Frailty Models}
\usage{tau(x)}
\arguments{
  \item{x}{A parametric frailty model, object of class \code{\link{parfm}}.}
}
%% \details{}
%% \value{}
\references{
  Marco Munda, Federico Rotolo, Catherine Legrand (2012).
  parfm: Parametric Frailty Models in R.
  Journal of Statistical Software, 51(11), 1-20. 
  URL http://www.jstatsoft.org/v51/i11/.

  Duchateau L., Janssen P. (2008) The frailty model. Springer.

  Wienke A. (2010) Frailty Models in Survival Analysis. Chapman & Hall/CRC
biostatistics series. Taylor and Francis.
}
\author{Federico Rotolo <federico.rotolo@stat.unipd.it> and Marco Munda
<marco.munda@uclouvain.be>}
%% \note{}

\seealso{
  \code{\link{parfm}}
}

\examples{
library(parfm)

#------Kidney dataset------
data(kidney) 
 # type 'help(kidney)' for a description of the data set
kidney$sex <- kidney$sex - 1

mod <- parfm(Surv(time,status) ~ sex + age, cluster="id",
             data=kidney, dist="exponential", frailty="gamma")
tau(mod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Kendall}
\keyword{tau}
\keyword{survival}
\keyword{frailty}
\keyword{parametric}
\keyword{marginal likelihood}
\keyword{shared}
