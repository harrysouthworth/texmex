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
dgpd(x, sigma, xi, u = 0, log.d = FALSE)
pgpd(q, sigma, xi, u = 0, lower.tail = TRUE, log.p = FALSE)
qgpd(p, sigma, xi, u = 0, lower.tail = TRUE, log.p = FALSE)
qgpd2(p, sigma, xi, u = 0, la = 1)
rgpd(n, sigma, xi, u = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p}{ Value or quantile. }
  \item{n}{ Number of random numbers to simulate.}
  \item{sigma}{ Scale parameter. }
  \item{xi}{ Shape parameter. }
  \item{u}{ Threshold }
  \item{la}{Rate at which threshold is exceeded.}
  \item{log.d, log.p}{Whether or not to work on the log scale.}
  \item{lower.tail}{Whether to return the lower tail.}
}
\details{
  The functions were originally based on code due to Stuart Coles and which
  appears in the \code{ismev} package. The functions have been vectorized and
  tested.

  \code{qgpd2} is used internally by \code{rl.plot} and is unlikely to be
  needed for any other reason.

  Random number generation is done by inversion of the distribution function.
  Code written by Harry Southworth.
}


\author{ Janet E Heffernan, Harry Southworth }
\examples{
  x <- rgpd(1000, sigma=1, xi=.5)
  hist(x)
  x <- rgpd(1000, sigma=exp(rnorm(1000, 1, .25)), xi=rnorm(1000, .5, .2))
  hist(x)
  plot(pgpd(x, sigma=1, xi=.5))
}
\keyword{ models }

