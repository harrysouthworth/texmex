\name{endPoint}
\alias{endPoint}
\alias{endPoint.evmOpt}
\alias{endPoint.evmSim}
\alias{endPoint.evmBoot}
\title{ Calculate upper end point for fitted extreme value model}
\description{
  Calculate upper end point for fitted extreme value model
}
\usage{
endPoint(y, verbose=TRUE, .unique=TRUE, ...)

\method{endPoint}{evmOpt}(y,verbose=TRUE, .unique=TRUE, ...)
\method{endPoint}{evmSim}(y,verbose=TRUE, .unique=TRUE, ...)

}
\arguments{
  \item{y}{ Object of class \code{evmOpt} or \code{evmSim}.}
  \item{verbose}{Whether to print output.}
  \item{.unique}{Whether or not to use only unique values of \code{y}.}
  \item{...}{further arguments to be passed to the \code{\link{signif}} function.}
}
\value{
  In cases where the fitted shape parameter is negative, the fitted finite upper endpoint of the GPD.
}

\author{ Janet E. Heffernan}

