\name{parfm}
\alias{parfm}
\alias{coef.parfm}
\alias{logLik.parfm}
\alias{plot.parfm}
\alias{print.parfm}
\title{Parametric Frailty Models}
\description{Fits parametric Frailty Models by maximizing the marginal likelihood.}
\encoding{UTF-8}

\usage{
parfm(formula, cluster = NULL, strata = NULL, data,
      inip = NULL, iniFpar = NULL,
      dist = c("weibull", "inweibull", "frechet", "exponential", 
               "gompertz", "loglogistic", "lognormal", "logskewnormal"),
      frailty   = c("none", "gamma", "ingau", "possta",
                    "lognormal", "loglogistic"),
      method = "nlminb", 
      maxit = 500, Fparscale = 1, showtime = FALSE, correct = 0)
}

\arguments{
  \item{formula}{A \code{\link{formula}} object, with the response on the left of a ~
operator,
and the terms on the right.
    The response must be a survival object as returned by the \code{\link{Surv}()}
function.
    The status indicator 'event' in the Surv object must be 0=alive, 1=dead.}
  \item{cluster}{The name of a cluster variable in data.}
  \item{strata}{The name of a strata variable in data.}
  \item{data}{A \code{\link{data.frame}} in which to interpret the variables named in
the formula.}
  \item{inip}{The vector of initial values.
    First components are for the baseline hazard parameters according to
the order
given in 'details';
    Other components are for the regression parameters according to the
order given
in 'formula'.}
  \item{iniFpar}{The initial value of the frailty parameter.}
  \item{dist}{
    The baseline hazard distribution.
    One of \code{weibull}, \code{inweibull}, \code{frechet}, \code{exponential},
           \code{gompertz}, \code{lognormal}, \code{loglogistic}, or \code{logskewnormal}.}
  \item{frailty}{
    The Frailty distribution.
    One of: \code{none}, \code{gamma}, \code{ingau}, \code{possta} or \code{lognormal}.}
  \item{method}{The optimisation method from the function \code{\link{optimx}()}.}
  \item{maxit}{Maximum number of iterations (see \code{\link{optimx}()}).}
  \item{Fparscale}{the scaling value for the frailty parameter in \code{\link{optimx}()}. 
    Optimisation is performed on Fpar/Fparscale.}
  \item{showtime}{Show the execution time?}
  \item{correct}{A correction factor that does not change the marginal
    log-likelihood except for an additive constant given by
    #clusters * correct * log(10). It may be useful in order to get finite
    log-likelihood values in case of many events per cluster
    with Positive Stable frailties. Note that the value of the
    log-likelihood in the output is the re-adjusted value.}
}

\details{

\strong{Baseline hazards}

The Weibull hazard (\code{dist="weibull"}) is
\deqn{h(t; \rho, \lambda) = \rho \lambda t^{\rho - 1}}{
  h(t; rho, lambda) = rho * lambda * t^(rho - 1)}
with \eqn{\rho, \lambda > 0}{rho, lambda > 0}.

The inverse Weibull (or Fréchet) hazard (\code{dist="inweibull"} or \code{dist="frechet"}) is
\deqn{h(t; \rho, \lambda) = \frac{\rho \lambda t^{-\rho - 1}}{\exp(\lambda t^{-\rho}) - 1}}{
  h(t; rho, lambda) = rho * lambda * t^(-rho - 1) / (exp(lambda t^-rho) - 1)}
with \eqn{\rho, \lambda > 0}{rho, lambda > 0}.

The exponential hazard (\code{dist="exponential"}) is
\deqn{h(t; \lambda) = \lambda}{
  h(t; lambda) = lambda}
with \eqn{\lambda > 0}{lambda > 0}.

The Gompertz hazard (\code{dist="gompertz"}) is
\deqn{h(t; \gamma, \lambda) = \lambda \exp ( \gamma t )}{
  h(t; gamma, lambda) = lambda * exp(gamma * t)}
with \eqn{\gamma, \lambda > 0}{ gamma, lambda > 0}.

The lognormal hazard (\code{dist="lognormal"}) is
\deqn{h(t; \mu, \sigma) = \frac{\phi \left( \frac{\log ( t ) - \mu}{\sigma}\right)}{\sigma t \left[ 1 - \Phi \left( \frac{\log ( t ) - \mu}{\sigma}\right)\right]}}{
  h(t; mu, sigma) = phi ( (log(t) - mu) / sigma ) / {sigma * t * [1 - Phi( (log(t)- mu) / sigma ) )]}}
with \eqn{\mu \in   {R}, \sigma > 0}{mu in R, sigma > 0}, and where
\eqn{\phi}{phi} and \eqn{\Phi}{Phi} are the density and distribution functions of a standard Normal random variable.

The loglogistic hazard (\code{dist="loglogistic"}) is
\deqn{h(t; \alpha, \kappa) = \frac{\exp ( \alpha ) \kappa t^{\kappa - 1}}{1 + \exp ( \alpha ) t^{\kappa}}}{
  h(t; alpha, kappa) = [exp(alpha) * kappa * t^(kappa - 1)] / [1 + exp(alpha) * t^(kappa)]}
with \eqn{\alpha \in   {R}, \kappa > 0}{alpha in R, kappa > 0}.

The log-skew-normal hazard (\code{dist="logskewnormal"}) is 
    obtained as the ratio between the density and the cumulative distribution function
    of a log-\link[sn:dsn]{skew normal} random variable (Azzalini, 1985),
    which has density
\deqn{f(t; \xi, \omega, \alpha) = \frac{2}{\omega t}
    \phi\left(\frac{\log(t) - \xi}{\omega}\right)
    \Phi\left(\alpha\frac{\log(t)-\xi}{\omega}\right)}
with \eqn{\xi \in {R}, \omega > 0, \alpha \in {R}}{xi in R, omega > 0, alpha in R}, and where
\eqn{\phi}{phi} and \eqn{\Phi}{Phi} are the density and distribution functions of a standard Normal random variable.
Of note, if \eqn{\alpha=0}{alpha=0} then the log-skew-normal boils down 
    to the log-normal distribution, with \eqn{\mu=\xi}{mu=xi} and \eqn{\sigma=\omega}{sigma=omega}.


\strong{Frailty distributions}

The gamma frailty distribution (\code{frailty="gamma"}) is
\deqn{f ( u ) = \frac{\theta^{-\frac{1}{\theta}} u^{\frac{1}{\theta} - 1} \exp \left( - u / \theta \right)}	{\Gamma ( 1 / \theta )}}{
f(u) = [theta^(-1 / theta) * u^(1 / theta - 1) * exp(-u / theta)] / [Gamma(1 / theta)]}
with \eqn{\theta > 0}{theta > 0} and where \eqn{\Gamma ( \cdot )}{Gamma(.)} is the gamma function.

The inverse Gaussian frailty distribution (\code{frailty="ingau"}) is
\deqn{f(u) = \frac1{\sqrt{2 \pi \theta}} u^{- \frac32} \exp \left( - \frac{(u-1)^2}{2 \theta u}  \right)}{ 
  f(u) = 1 / sqrt(2 * pi * theta) * u^(-3 / 2) * exp[- (u - 1)^2 / (2 * theta * u)]}
with \eqn{\theta > 0}{theta > 0}.

The positive stable frailty distribution (\code{frailty="poosta"}) is
\deqn{f(u) = f(u) = - \frac1{\pi u} \sum_{k=1}^{\infty} \frac{\Gamma ( k (1 - \nu ) + 1 )}{k!} \left( - u^{ \nu - 1} \right)^{k} \sin ( ( 1 - \nu ) k \pi )}{
f(u) = - 1 / (pi * u) * sum_(k=1, 2, ..., infiny) [Gamma(k * (1 - nu) + 1) / k! * (-u^(nu - 1))^k * sin((1 - nu) * k * pi)]}
with \eqn{0 < \nu < 1}{0 < nu < 1}.

The lognormal frailty distribution (\code{frailty="lognormal"}) is
\deqn{f(u) = \frac1{\sqrt{2 \pi \theta}} u^{-1} \exp \left( - \frac{\log(u)^2}{2 \theta}  \right)}{ 
  f(u) = 1 / sqrt(2 * pi * theta) * u^(-1) * exp[- (log u)^2 / (2 * theta)]}
with \eqn{\theta > 0}{theta > 0}.
As the Laplace tranform of the lognormal frailties does not exist in closed form,
  the saddlepoint approximation is used (Goutis and Casella, 1999).

}

\value{An object of class \code{parfm}.}

\references{
  Munda M, Rotolo F, Legrand C (2012).
  parfm: Parametric Frailty Models in R.
  Journal of Statistical Software, 51(11), 1-20. 
  DOI <\doi{10.18637/jss.v051.i11}>

  Duchateau L, Janssen P (2008). The frailty model. Springer.

  Wienke A (2010).
  Frailty Models in Survival Analysis. 
  Chapman & Hall/CRC biostatistics series. Taylor and Francis.

  Goutis C, Casella G (1999).
  Explaining the Saddlepoint Approximation.
  The American Statistician 53(3), 216-224.
  DOI <\doi{10.1080/00031305.1999.10474463}>
    
 Azzalini A (1985).
 A class of distributions which includes the normal ones.
 Scandinavian Journal of Statistics, 12(2):171-178.
}
\author{
  \packageAuthor{parfm}
}


\seealso{
  \code{\link{select.parfm}},
  \code{\link{ci.parfm}},
  \code{\link{predict.parfm}}
}

\examples{
#------Kidney dataset------
data(kidney) 
 # type 'help(kidney)' for a description of the data set
kidney$sex <- kidney$sex - 1

parfm(Surv(time,status) ~ sex + age, cluster = "id",
      data = kidney, dist = "exponential", frailty = "gamma")
      
      
\donttest{
parfm(Surv(time,status) ~ sex + age, cluster = "id",
      data = kidney, dist = "exponential", frailty = "lognormal")

parfm(Surv(time,status) ~ sex + age, cluster = "id",
      data = kidney, dist = "weibull", frailty = "ingau")

parfm(Surv(time,status) ~ sex + age, cluster = "id",
      data = kidney, dist="gompertz", frailty="possta", method="CG")


parfm(Surv(time,status) ~ sex + age, cluster = "id",
      data = kidney, dist="logskewnormal", frailty="possta", method = 'BFGS')
#--------------------------

#------Asthma dataset------
data(asthma)
head(asthma)
  # type 'help(asthma)' for a description of the data set

asthma$time <- asthma$End - asthma$Begin
parfm(Surv(time, Status) ~ Drug, cluster = "Patid", data = asthma,
      dist = "weibull", frailty = "gamma")
      
parfm(Surv(time, Status) ~ Drug, cluster = "Patid", data = asthma,
      dist = "weibull", frailty = "lognormal")

parfm(Surv(Begin, End, Status) ~ Drug, cluster = "Patid", 
      data = asthma[asthma$Fevent == 0, ],
      dist = "weibull", frailty = "lognormal", method = "Nelder-Mead")
#--------------------------
}
}

\keyword{survival}
\keyword{frailty}
\keyword{parametric}
\keyword{marginal likelihood}
\keyword{shared}
