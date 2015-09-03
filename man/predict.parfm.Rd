\name{predict.parfm}
\alias{predict.parfm}
\alias{plot.predict.parfm}
\alias{print.predict.parfm}
\title{Predictions of frailty values for Parametric Frailty Models}

\description{The function \code{predict.parfm()} computes predictions of frailty values for objects of class \code{\link{parfm}}.}

\usage{
\method{predict}{parfm}(object, ...)
}

\arguments{
  \item{object}{A parametric frailty model, object of class \code{\link{parfm}}.}
  \item{...}{see \code{\link{predict}()}}
}

%%\details{}
\value{An object of class \code{predict.parfm}.}
\references{
  Glidden D., Vittinghoff E. (2004). 
  Modelling Clustered Survival Data From Multicentre Clinical Trials.   
  Statistics in medicine, 23(3), 369--388.

  Marco Munda, Federico Rotolo, Catherine Legrand (2012).
  parfm: Parametric Frailty Models in R.
  Journal of Statistical Software, 51(11), 1-20. 
  URL http://www.jstatsoft.org/v51/i11/.
  
%  McGilchrist C.A., Aisbett C.W. (1991). 
%  Regression with Frailty in Survival Analysis. Biometrics, pp. 461--466.
}

\author{Federico Rotolo <federico.rotolo@stat.unipd.it> and Marco Munda <marco.munda@arlenda.be>}

%%\note{}

\seealso{\code{\link{parfm}}}

\examples{
library(parfm)
data(kidney)
kidney$sex <- kidney$sex - 1

model <- parfm(Surv(time,status) ~ sex + age, 
               cluster="id", data=kidney,
	             dist="exponential", frailty="gamma")
u <- predict(model)
u


# Predictions from semi-parametric Gamma frailty model
# via coxph() function
model.coxph <- coxph(Surv(time,status) ~ sex + age + 
                     frailty(id, frailty="gamma", eps=1e-11), 
                     outer.max=15, data=kidney)
u.coxph <- exp(model.coxph$frail)


# Plot of predictions from both models
par(mfrow=c(1,2))
ylim <- c(0, max(c(u, u.coxph)))
plot(u, sort="i",
     main=paste("Parametric", 
                "Gamma frailty model",
                "with Exponential baseline", 
                sep="\n"),
     ylim=ylim)
                
names(u.coxph) <- kidney[seq(2,76, 2), "id"]
class(u.coxph) <- "predict.parfm"
attr(u.coxph, "clustname") <- "id"
plot(u.coxph, sort="i",
    main=paste("Semi-parametric",
               "Gamma frailty model", sep="\n"),
    ylim=ylim)
}
\keyword{frailty}
\keyword{prediction}
\keyword{frailty prediction}
