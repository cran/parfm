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
  Kendall, MG. (1938)
  A new measure of rank correlation.
  Biometrika 30(1/2), 81-93.
  DOI:   DOI \href{http://dx.doi.org/10.2307/2332226}{10.2307/2332226}
}
\author{
  \packageAuthor{parfm}
}

\seealso{
  \code{\link{parfm}}
}

\examples{
library(parfm)

data(kidney) 
 # type 'help(kidney)' for a description of the data set
kidney$sex <- kidney$sex - 1

mod <- parfm(Surv(time,status) ~ sex + age, cluster = "id",
             data = kidney, dist = "exponential", frailty = "gamma")
tau(mod)
}


\keyword{Kendall}
\keyword{tau}
\keyword{survival}
\keyword{frailty}
\keyword{parametric}
\keyword{marginal likelihood}
\keyword{shared}
