\name{select.parfm}
\alias{select.parfm}
\alias{plot.select.parfm}
\alias{print.select.parfm}
\title{AIC and BIC values of several Parametric Frailty Models}

\description{The function \code{select.parfm()} computes the AIC and BIC values
  of parametric frailty models with different baseline hazards and different frailty distributions.}

\usage{
select.parfm(formula, cluster=NULL, strata=NULL, data, inip=NULL, iniFpar=NULL,
             dist=c("exponential", "weibull", "inweibull", "frechet", "gompertz", 
                    "loglogistic", "lognormal", "logskewnormal"),
             frailty=c("none", "gamma", "ingau", "possta", "lognormal"),
             method="BFGS", maxit=500, Fparscale=1, correct=0)
}

\arguments{
  \item{formula}{A \code{\link{formula}} object, with the response on the left of a ~ operator, and the terms on the right.
    The response must be a survival object as returned by the \code{\link{Surv}} function.}
  \item{cluster}{The name of a cluster variable in data.}
  \item{strata}{The name of a strata variable in data.}
  \item{data}{A \code{\link{data.frame}} in which to interpret the variables named in the formula.}
  \item{inip}{The vector of initial values.
    First components are for the baseline hazard parameters according to
the order
given in 'details';
    Other components are for the regression parameters according to the
order given
in 'formula'.}
  \item{iniFpar}{The initial value of the frailty parameter.}
  \item{dist}{The vector of baseline hazards' names.
    It can include any of \code{weibull}, \code{inweibull} (alias \code{frechet}),
        \code{exponential}, \code{gompertz}, \code{loglogistic} or \code{lognormal}.}
  \item{frailty}{The vector of  frailty distributions' names.
    It can include any of: \code{none}, \code{gamma}, \code{ingau}, \code{possta} or \code{lognormal}.}
  \item{method}{The optimisation method from the function \code{\link{optim}()}.}
  \item{maxit}{Maximum number of iterations (see \code{\link{optim}()}).}
  \item{Fparscale}{the scaling value for the frailty parameter in \code{\link{optim}()}. 
    Optimisation is performed on Fpar/Fparscale.}
  \item{correct}{A correction factor that does not change the marginal
    log-likelihood except for an additive constant given by
    #clusters * correct * log(10). It may be useful in order to get finite
    log-likelihood values in case of many events per cluster
    with Positive Stable frailties. Note that the value of the
    log-likelihood in the output is the re-adjusted value.}
}

%%\details{}
\value{An object of class \code{select.parfm}.}

\author{
  \packageAuthor{parfm}
}

%%\note{}
\references{
  Munda M, Rotolo F, Legrand C (2012).
  parfm: Parametric Frailty Models in R.
  Journal of Statistical Software, 51(11), 1-20. 
  DOI <\doi{10.18637/jss.v051.i11}>
}


\seealso{
  \code{\link{parfm}},
  \code{\link{ci.parfm}},
  \code{\link{predict.parfm}}
}

\examples{
\donttest{
data(kidney)
kidney$sex <- kidney$sex - 1

models <- select.parfm(Surv(time,status) ~ sex + age, 
                       dist = c("exponential", 
                                "weibull",
                                "inweibull",
                                "loglogistic", 
                                "lognormal", 
                                "logskewnormal"),
                       frailty = c("gamma", 
                                   "ingau", 
                                   "possta", 
                                   "lognormal"),
                       cluster = "id", data = kidney)
models
plot(models)
}
}
%%\keyword{}
