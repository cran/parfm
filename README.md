[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/parfm)](https://cran.r-project.org/package=parfm)

# `parfm`, Parametric Frailty Models in `R`
Federico Rotolo and Marco Munda

## Description
Fits Parametric Frailty Models by maximum marginal likelihood.
             Possible baseline hazards:
                 exponential, Weibull, inverse Weibull (Fréchet),
                 Gompertz, lognormal, log-skew-normal, and loglogistic.
             Possible Frailty distributions:
                gamma, positive stable, inverse Gaussian and lognormal.
                
## Details
Frailty models are survival models for clustered or overdispersed time-to-event data.
They consist in proportional hazards Cox's models with the addition of a random effect,
  accounting for different risk levels.

When the form of the baseline hazard is somehow known in advance, the parametric estimation approach can be used advantageously.
The `parfm` package provides a wide range of parametric frailty models in `R`.
The following baseline hazard families are implemented

* exponential,

* Weibull,

* inverse Weibull (Fréchet),

* Gompertz,

* lognormal,

* log-skew-normal,

* loglogistic,

together with the frailty distributions 

* gamma,

* positive stable,

* inverse Gaussian, and 

* lognormal.

Parameter estimation is done by maximising the marginal log-likelihood,
  with right-censored and possibly left-truncated data.


## Parametrisations
### Baseline hazards
The **exponential** hazard is
  $$h(t; \lambda) = \lambda,$$
with $\lambda > 0$.

The **Weibull** hazard is
  $$h(t; \rho, \lambda) = \rho \lambda t^{\rho-1},$$
with $\rho,\lambda > 0$.

The **inverse Weibull** (or **Fréchet**) hazard is
  $$h(t; \rho, \lambda) = \frac{\rho \lambda t^{-\rho - 1}}{\exp(\lambda t^{-\rho}) - 1}$$
with $\rho, \lambda > 0$.

  $$h(t; \rho, \lambda) = \rho \lambda t^{\rho-1},$$
with $\rho,\lambda > 0$.

The **Gompertz** hazard is
  $$h(t; \gamma, \lambda) = \lambda e^{\gamma t},$$
with $\gamma,\lambda > 0$.

The **lognormal** hazard is
  $$h(t; \mu, \sigma) = 
    { \phi([log t -\mu]/\sigma)} / { \sigma t [1-\Phi([log t -\mu]/\sigma)]},$$
with $\mu\in\mathbb R$, $\sigma > 0$ and $\phi(\cdot)$ and $\Phi(\cdot)$
the density and distribution functions of a standard Normal.


The **log-skew-normal** hazard is obtained as the ratio between the density
  and the cumulative distribution function
  of a log-skew normal random variable (Azzalini, 1985),
  which has density
$$f(t; \xi, \omega, \alpha) = \frac{2}{\omega t}
    \phi\left(\frac{\log(t) - \xi}{\omega}\right)
    \Phi\left(\alpha\frac{\log(t)-\xi}{\omega}\right)$$
with $\xi \in {R}, \omega > 0, \alpha \in {R}$, and where
$\phi(\cdot)$ and $\Phi(\cdot)$ are the density and distribution functions of a standard Normal random variable.
Of note, if $alpha=0$ then the log-skew-normal boils down 
to the log-normal distribution, with $\mu=\xi$ and $\sigma=\omega$.
    
The **loglogistic** hazard is
  $$h(t; \alpha, \kappa) = 
    {exp(\alpha) \kappa t^{\kappa-1} } / {
      1 + exp(\alpha) t^{\kappa}},$$
with $\alpha\in\mathbb R$ and $\kappa>0$.

### Frailty distributions
The **gamma** frailty distribution is
$$f ( u ) = \frac{\theta^{-\frac{1}{\theta}} u^{\frac{1}{\theta} - 1} \exp \left( - u / \theta \right)}	{\Gamma ( 1 / \theta )}$$
with $\theta > 0$ and where $\Gamma(\cdot)$ is the gamma function.

The **inverse Gaussian** frailty distribution is
$$f(u) = \frac1{\sqrt{2 \pi \theta}} u^{- \frac32} \exp \left( - \frac{(u-1)^2}{2 \theta u}  \right)$$
with $\theta > 0$.

The **positive stable** frailty distribution is
$$f(u) = f(u) = - \frac1{\pi u} \sum_{k=1}^{\infty} \frac{\Gamma ( k (1 - \nu ) + 1 )}{k!} \left( - u^{ \nu - 1} \right)^{k} \sin ( ( 1 - \nu ) k \pi )$$
with $0 < \nu < 1$.

The **lognormal** frailty distribution is
$$f(u) = \frac1{\sqrt{2 \pi \theta}} u^{-1} \exp \left( - \frac{\log(u)^2}{2 \theta}  \right)$$
with $\theta > 0$.
As the Laplace tranform of the lognormal frailties does not exist in closed form,
  the saddlepoint approximation is used (Goutis and Casella, 1999).

---
## References
    
Azzalini A (1985).
 A class of distributions which includes the normal ones.
 *Scandinavian Journal of Statistics*, **12**(2):171-178.
 URL [http://www.jstor.org/stable/4615982]


Cox DR (1972).
  Regression models and life-tables.
  *Journal of the Royal Statistical Society. Series B (Methodological)*, **34**:187–220.

Duchateau L, Janssen P (2008). *The frailty model*. Springer.

Goutis C, Casella G (1999).
  Explaining the Saddlepoint Approximation.
  *The American Statistician*, **53**(3):216-224.
[10.1080/00031305.1999.10474463](http://dx.doi.org/10.1080/00031305.1999.10474463).

Munda M, Rotolo F, Legrand C (2012).
  parfm: Parametric Frailty Models in R.
  *Journal of Statistical Software*, **51**(11):1-20. 
  DOI: [10.18637/jss.v051.i11](http://dx.doi.org/10.18637/jss.v051.i11)

Wienke A (2010).
  *Frailty Models in Survival Analysis*. 
  Chapman & Hall/CRC biostatistics series. Taylor and Francis.