\name{gpd}
\alias{gpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit a generlized Pareto distribution by maximum
    (penalized) likelihood }
\description{
  Fit a generlized Pareto distribution by maximum
    (penalized) likelihood, possibly with explanatory variables.
}
\usage{
gpd(y, data, th, phi = ~1, xi = ~1, penalty = "gaussian", xi.var = 10^4, gaussianParameters = NULL, maxit = 10000, trace = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Either a numeric vector or the name of a variable in \code{data}. }
  \item{data}{ A data object containing \code{y} and any covariates. }
  \item{th}{ The threshold for \code{y} }
  \item{phi}{ Formula for the log of the scale parameter. Defaults to
      \code{phi = ~ 1} - i.e. no covariates. }
  \item{xi}{ Formula for the shape parameter. Defaults to
      \code{phi = ~ 1} - i.e. no covariates. }
  \item{penalty}{ How to penalize the likelhood. Currently, either
     ``none'' or ``gaussian'' are the only allowed arguments. }
  \item{xi.var}{ XXX CHECK IF THIS IS USED }
  \item{gaussianParameters}{  A list with two components. The first should be 
                a vector of means, the second should be a covariance matrix.
                These represent the parameters in the Gaussian prior distribution.
                If not supplied: all default prior means are zero;
                all default prior variances are $10^4$;all covariances
                are zero.  }
  \item{maxit}{ The number of iterations allowed in \code{optim} }
  \item{trace}{ Whether or not to print progress to screen. }
}
\details{
  Working with the log of the scale parameter was found to improve the
  stability of computations, and it makes a quadratic penalty more appropriate.
  A quadratic penalty can be thought of as a Gaussian prior distribution, whence
  the terminology of the function.
  
  Some of the internal code has been taken from the \code{gpd.fit} function
  in the \code{ismev} package and is due to Janet E Heffernan.
}
\value{
  A list containing the model.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Janet E Heffernan, Harry Southworth }

\examples{
  x <- rnorm(1000)
  mod <- gpd(x, th=quantile(x, .7))
  mod
}
\keyword{ models }

