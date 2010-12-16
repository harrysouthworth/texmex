\name{copula}
\alias{copula}
\alias{plot.copula}
\title{ Get the copula of a matrix of variables }
\description{
Returns the copula of several random variables.
}
\usage{
copula(x, na.last = NA)
\method{print}{copula}(object, ...)
\method{show}{copula}(object, ...)
\method{summary}{copula}(object, ...)
\method{plot}{copula}(object, jitter=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix or data.frame containing numeric variables. }
  \item{na.last}{ How to treat missing values. See \code{rank} for details. }
  \item{jitter}{In the call to \code{plot.copula}, if \code{jitter=TRUE},
                the values are jittered before plotting. Defaults
                to \code{jitter=FALSE}.}
  \item{object}{Object of class \code{copula} as returned by function \code{copula}.}
  \item{...}{Further arguments to be passed to methods.}                
}
\details{
  The result is obtained by applying \code{\link{edf}} to each
  column of \code{x} in turn.
}
\value{
  A matrix with the same dimensions as \code{x}, each column
  of which contains the quantiles of each column of \code{x}.
}
\author{ Harry Southworth }
\seealso{ \code{\link{edf}} }
\examples{
  D <- liver[liver$dose == "D",]
  Dco <- copula(D)
  plot(Dco)
}

\keyword{ multivariate }

