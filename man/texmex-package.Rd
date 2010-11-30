\name{texmex-package}
\alias{texmex-package}
\alias{texmex}
\docType{package}
\title{
Conditional multivariate extreme value modelling
}
\description{
Modelling extreme values using generalized Pareto distributions
and the Heffernan-Tawn appraoch to conditional multivariate extremes
}
\details{
\tabular{ll}{
Package: \tab texmex\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2010-06-22\cr
License: \tab GPL (>=2) | BSD\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
}
Fit generalized Pareto distributions by (penalized) maximum likelihood
(\code{gpd}), perform Bayesian inference via Markov chain Monte Carlo
(\code{bgpd}), perform conditional multivariate extreme value modelling
following the approach of Heffernan and Tawn (\code{mexDependece},
\code{bootmex}, \code{mexDependence}).
}
\author{
Harry Southworth

Maintainer: Harry Southworth <harry.southworth@astrazeneca.com>

}
\references{

S. Coles, An Introduction to Statistical Modeling of Extreme Values,
    Springer, 2001
    
J. E. Heffernan and J. A. Tawn, A conditional approach for multivariate
    extreme values, Journal of the Royal Statistical Society B, 66, 497 -- 546, 
    2004
}

\keyword{ package }

\examples{

alt <- log(liver$ALT_M[liver$dose == "D"])

# Do plots to aid threshold selection
par(mfrow=c(2, 2))
mrlPlot(alt)
gpdRangeFit(alt, umin=quantile(alt, .25), umax=quantile(alt, .9))

mod <- gpd(alt, th=quantile(alt, .7))
summary(mod)
bmod <- bgpd(alt, th=quantile(alt, .7))
par(mfrow=c(3, 2))
plot(bmod)

D <- liver[liver$dose == "D", c("ALT_M", "AST_M", "ALP_M", "TBL_M")]
gpdD <- migpd(D, qu=.7) # Use data over 70th percentile
summary(gpdD)
mexDependence(gpdD)
bootD <- bootmex(gpdD, which="ALT_M", gqu=.7)
bootD
predD <- predict.mex(gpdD, bootD)
predD
par(mfrow=c(2, 2))
plot(predD)
}
