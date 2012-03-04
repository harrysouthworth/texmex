\name{copula}
\alias{copula}
\alias{plot.copula}
\title{ Calculate the copula of a matrix of variables }
\description{
Returns the copula of several random variables.
}
\usage{
copula(x, na.last = NA)
\method{plot}{copula}(x, jitter. = FALSE, jitter.factor=1, ...)
}
\arguments{
  \item{x}{ A matrix or data.frame containing numeric variables. }
  \item{na.last}{ How to treat missing values. See \code{rank} for details. }
  \item{jitter.}{In the call to \code{plot.copula}, if \code{jitter=TRUE},
                the values are jittered before plotting. Defaults
                to \code{jitter. = FALSE}.}
  \item{jitter.factor}{How much jittering to use. Defaults to \code{jitter.factor = 1.}}
  \item{...}{Further arguments to be passed to plot method.}                
}
\details{
  The result is obtained by applying \code{\link{edf}} to each
  column of \code{x} in turn.
  
  Print, plot, show and summary methods are available for the copula class.
}
\value{
  A matrix with the same dimensions as \code{x}, each column
  of which contains the quantiles of each column of \code{x}. This object uis of class \code{copula}.
}
\author{ Harry Southworth }
\seealso{ \code{\link{edf}} }
\examples{
  D <- liver[liver$dose == "D",]
  Dco <- copula(D)
  plot(Dco)
}

\keyword{ multivariate }

