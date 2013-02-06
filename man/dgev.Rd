\name{dgev}
\alias{dgev}
\alias{pgev}
\alias{qgev}
\alias{rgev}
\title{ Density, cumulative density, quantiles and random number generation
        for the generalized extreme value distribution }
\description{
  Density, cumulative density, quantiles and random number generation
        for the generalized Pareto distribution
}
\usage{
dev(x, sigma, xi, u = 0, log.d = FALSE)
pev(q, sigma, xi, u = 0, lower.tail = TRUE, log.p = FALSE)
qev(p, sigma, xi, u = 0, lower.tail = TRUE, log.p = FALSE)
rev(n, sigma, xi, u = 0)
}

\arguments{
  \item{x, q, p}{ Value, quantile or probability respectively. }
  \item{n}{ Number of random numbers to simulate.}
  \item{sigma}{ Scale parameter. }
  \item{xi}{ Shape parameter. }
  \item{u}{ Threshold }
  \item{log.d, log.p}{Whether or not to work on the log scale.}
  \item{lower.tail}{Whether to return the lower tail.}
}
\details{
  Random number generation is done by inversion of the distribution function.
}


\author{Harry Southworth}
\examples{
  x <- rgev(1000, sigma=1, xi=.5)
  hist(x)
  x <- rgev(1000, sigma=exp(rnorm(1000, 1, .25)), xi=rnorm(1000, .5, .2))
  hist(x)
  plot(pgev(x, sigma=1, xi=.5))
}
\keyword{ models }

