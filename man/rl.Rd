\name{rl}
\alias{rl}
\alias{rl.gpd}
\alias{rl.bgpd}
\alias{rl.bootgpd}
\docType{data}
\title{
Return levels
}
\description{
Computation of return levels and confidence intervals for generalized Pareto distributions.
}
\usage{rl(object, M=1000, newdata)
\method{rl}{gpd}(object, M=1000, newdata)
\method{rl}{bgpd}(object, M=1000, newdata)
\method{rl}{bootgpd}(object, M=1000, newdata)
}
\arguments{

\item{object}{An object of class \code{gpd}, \code{bgpd} or \code{bootgpd}.}
\item{M}{The M-observation return level is computed by the function. Defaults to \code{M = 1000}.}
\item{newdata}{Data from which to calculate the return level. If not provided, the original data
               used to fit the model is used.}
}

\details{
The M-observation return level is defined as the value that is expected to be exceeded only
once every M observations. Thus, it is an estimate of a high quantile of the fitted distribution.
In models fit by the \code{gbm} family of functions, only a fraction of the data is actually
included in the model. This consideration is taken into account by \code{rl}.
}
\examples{
mod <- gpd(rain, qu=.8)
rl(mod, M=100*365) # 100-year return level
}
\keyword{datasets}
