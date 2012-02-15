\name{ci.parfm}
\alias{ci.parfm}
\title{Confidence Intervals for Regression Parameter Estimates of Parametric Frailty Models}
\description{Computes confidence intervals of regression parameter estimates for objects of class \code{\link{parfm}}.}

\usage{
  ci.parfm(x, level=.05, digits=3)
}

\arguments{
  \item{x}{
    The fitted model, object of class \code{parfm}, obtained as result of the \code{\link{parfm}()} function.}
  \item{level}{
    The coverage probability of the confidence interval.}
  \item{digits}{
    The number of significant digits.}
}

%% \details{}
%%\value{}
%%\references{}

\author{Federico Rotolo <federico.rotolo@stat.unipd.it> and Marco Munda <marco.munda@uclouvain.be>}

%%\note{}

\seealso{
  \code{\link{parfm}},
  \code{\link{select.parfm}},
  \code{\link{predict.parfm}}
}

\examples{
  library(parfm)

  #------Kidney dataset------
  data("kidney")
  kidney$sex <- kidney$sex - 1

  pfm <- parfm(Surv(time,status) ~ sex + age, cluster="id", 
               data=kidney, dist="exponential", frailty="gamma")

  ci.parfm(pfm)
  #--------------------------
}

%%\keyword{}
