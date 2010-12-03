\name{MCS}
\alias{MCS}
\alias{bootMCS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multivariate conditional Spearman's rho }
\description{
  Compute multivariate conditional Spearman's rho over a
  range of quantiles.
}
\usage{
MCS(X, p = seq(0.1, 0.9, by = 0.1))
bootMCS(X, p = seq(0.1, 0.9, by = 0.1), R = 100, trace = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ A matrix of numeric variables. }
  \item{p}{ The quantiles at which to evaluate. }
  \item{R}{The number of bootstrap samples to run. Defaults to
          \code{R = 100}.}
  \item{trace}{How often to inform the user of progress. Defaults to
          \code{trace = 10}.}
}
\details{
  The method is described in detail by Schmid and Schmidt (2007).
  The main code was written by Yiannis Papastathopoulos, wrappers
  written by Harry Southworth.
  
  When the result of a call to \code{bootMCS} is plotted, simple
  quantile bootstrap confidence intervals are displayed.
}
\value{
  If it is a LIST, use
  \item{MCS }{The correlations.}
  \item{p }{The quantiles at which the correlations were evaluated at}
  \item{method}{The method used.}
  \item{replicates}{Bootstrap replicates if \code{bootMCS} was used.}
}
\references{ F. Schmid and R. Schmidt, Multivariate conditional versions
of Spearman's rho and related measures of tail dependence, Journal of
Multivariate Analysis, 98, 1123 -- 1140, 2007 }
\author{ Yiannis Papastathopoulos, Harry Southworth }

\seealso{ \code{\link{chi}} }
\examples{
D <- liver[liver$dose == "D",]
Dmcs <- bootMCS(D[, 5:6])
plot(Dmcs)

}
\keyword{ multivariate }

