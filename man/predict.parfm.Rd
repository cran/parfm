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

  McGilchrist C.A., Aisbett C.W. (1991). 
  Regression with Frailty in Survival Analysis. Biometrics, pp. 461--466.
  
%  Munda M., Rotolo F., Legrand C. (2012). parfm: Parametric Frailty Models in R.
% Submitted.
}

\author{Federico Rotolo <federico.rotolo@stat.unipd.it> and Marco Munda <marco.munda@uclouvain.be>}

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


par(mfrow=c(1,2))
plot(u, sort="i",
     main=paste("Parametric", "Gamma frailty model",
                "with Exponential baseline", sep="\n"))
                
# Predictions from the original paper
# McGilchrist and Aisbett (1991)
ukidney <- kidney[seq(2,76, 2), "frail"]
names(ukidney) <- kidney[seq(2,76, 2), "id"]
class(ukidney) <- "predict.parfm"
attr(ukidney, "clustname") <- "id"
plot(ukidney, sort="i",
     main=paste("Semi-parametric",
                "Lognormal frailty model", sep="\n"),
     hline=exp(.3821/2))
 # horizontal line at exp(theta/2),
 # i.e. the mean of Lognormal frailties
 # with log-frailties N(0, theta)

}
\keyword{frailty}
\keyword{prediction}
\keyword{frailty prediction}
