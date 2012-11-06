\name{anova.parfm}
\alias{anova.parfm}
\alias{anova.parfmlist}
\title{Analysis of Deviance for a parametric frailty model.}
\usage{
\method{anova}{parfm}(object, \dots)
}
\description{
  Compute an analysis of deviance table for one or more parametric frailty model fits.
}
\arguments{
  \item{object}{An object of class \code{parfm}}
  \item{\dots}{Further \code{parfm} objects}
}
\details{
  Specifying a single object gives a sequential analysis of deviance
  table for that fit.  That is, the reductions in the model log-likelihood
  as each term of the formula is added in turn are given in as
  the rows of a table, plus the log-likelihoods themselves.
  
  If more than one object is specified, the table has a row for the
  degrees of freedom and loglikelihood for each model. For all
  but the first model, the change in degrees of freedom and loglik
  is also given. (This only make statistical sense if the models are
  nested.)  It is conventional to list the models from smallest to
  largest, but this is up to the user.
  
  The table contains test statistics (and P values)
  comparing the reduction in loglik for each row.
}

\value{
  An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
}

\section{Warning}{
  The comparison between two or more models by \code{anova}
   will only be valid if they
  are fitted to the same dataset. This may be a problem if there are
  missing values.}
  
\seealso{
  \code{\link{parfm}}, \code{\link{anova}}.
}

\references{
  Marco Munda, Federico Rotolo, Catherine Legrand (2012).
  parfm: Parametric Frailty Models in R.
  Journal of Statistical Software, 51(11), 1-20. 
  URL http://www.jstatsoft.org/v51/i11/.
  }
  
\examples{
fit <- parfm(formula = Surv(time, status) ~ sex + age, cluster = "id",
             data = kidney, dist = "exponential", frailty = "gamma")
anova(fit)
fit2 <- parfm(formula = Surv(time, status) ~ sex, cluster = "id",
              data = kidney, dist = "exponential", frailty = "gamma")
anova(fit2,fit)
 }
\keyword{models}
\keyword{regression}
\keyword{survival}

