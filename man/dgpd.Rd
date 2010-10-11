\name{dgpd}
\alias{dgpd}
\alias{pgpd}
\alias{qgpd}
\alias{qgpd2}
\alias{rgpd}
\title{ Density, cumulative density, quantiles and random number generation
        for the generalized Pareto distribution }
\description{
  Density, cumulative density, quantiles and random number generation
        for the generalized Pareto distribution
}
\usage{
dgpd(x, sigma = 1, xi = 1, u = 0)
pgpd(q, sigma = 1, xi = 1, u = 0)
qgpd(p, sigma = 1, xi = 1, u = 0)
qgpd2(p, sigma = 1, xi = 1, u = 0, la = 1)
rgpd(n, sigma = 1, xi = 1, u = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p}{ Value or quantile. }
  \item{n}{ Number of random numbers to simulate.}
  \item{sigma}{ Scale parameter. }
  \item{xi}{ Shape parameter. }
  \item{u}{ Threshold }
}
\details{
  Most of the original code was written by Janet E Heffernan and 
  appears in the \code{ismev} package, reproduced here for convenience and
  to use more standard S/R function names.
  
  Random number generation is done by inversion of the distribution function.
  Code written by Harry Southworth.
}


\author{ Janet E Heffernan, Harry Southworth }
\examples{
  x <- rgpd(1000, sigma=1, xi=.5)
  hist(x)
  x <- rgpd(1000, sigma=exp(rnorm(1000, 1, .25)), xi=rnorm(1000, .5, .2))
  hist(x)
  plot(pgpd(x))
}
\keyword{ models }

