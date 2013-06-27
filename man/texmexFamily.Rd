\name{texmexFamily}
\alias{texmexFamily}
\alias{print.texmexFamily}
\alias{summary.texmexFamily}
\alias{gpd}
\alias{gev}

\title{Create families of distributions}

\description{Create families of distributions for use with extreme value modelling.}

\usage{
texmexFamily(name, log.lik, param, info = NULL, start = NULL, resid =
NULL, rl, delta, endpoint, density, rng, prob, quant)
\method{print}{texmexFamily}(x, ...)
\method{summary}{texmexFamily}(object, ...)
}

\arguments{
  \item{name}{The name of the distribution.}
  \item{log.lik}{The distribution's log-likelihood function.}
  \item{param}{The names of the parameters in the model.}
  \item{info}{Function to compute the information matrix. If not
    provided, the modelling functions will work with a numerical
    approximation.}
  \item{start}{Function to compute starting parameters for the
    model. If not provided, the modelling functions will try to
    guess.}
  \item{resid}{Function to compute residuals for the model.}
  \item{rl}{Function to compute return levels.}
  \item{delta}{Function to compute adjustments for covariance for return
    levels.}
  \item{endpoint}{Function to compute the upper or lower endpoint of the
    fitted distribution.}
  \item{density}{Function to compute the density.}
  \item{rng}{Function for random number generation.}
  \item{prob}{Function to compute cumulative probabilities.}
  \item{quant}{Function to compute quantiles.}
  \item{...}{Additional arguments to the print and summary methods.}
  \item{x, object}{An object of class 'texmexFamily'.}
}
\details{
  The \code{density}, \code{rng}, \code{prob} and \code{quant}
  functions can be simple wrappers for the usual d, r, p and q
  functions. They should take a matrix with number of columns
  equal to the number of parameters, and a fitted model object
  even if the model object is not used by the function.

  Examples of "texmexFamily" objects are \code{gpd} and \code{gev}.
  Take a look at those objects to see how the functions should
  be constructed.

  The functions are used by the modelling functions to create
  diagnostic plots, predictions, etc..
}
\value{
  A object of class "texmexFamily", which is essentially a list
  containing the input arguments. If
  \code{info}, \code{start}, \code{resid} are not provided, they
  default to \code{NULL}.
}
\author{Harry Southworth}
\note{
  The \code{gpd} and \code{gev} families are provided. The
  \code{\link{evm}} function defaults to using the \code{gpd}
  family.
}
\seealso{\code{\link{evm}}}
\keyword{ models }

