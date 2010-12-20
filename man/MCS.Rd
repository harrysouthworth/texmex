\name{MCS}
\alias{MCS}
\alias{bootMCS}
\alias{print.MCS}
\alias{show.MCS}
\alias{plot.MCS}
\alias{summary.MCS}
\alias{print.bootMCS}
\alias{show.bootMCS}
\alias{summary.bootMCS}
\alias{plot.bootMCS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multivariate conditional Spearman's rho }
\description{
  Compute multivariate conditional Spearman's rho over a
  range of quantiles.
}
\usage{
MCS(X, p = seq(0.1, 0.9, by = 0.1))
bootMCS(X, p = seq(0.1, 0.9, by = 0.1), R = 100, trace = 10)
\method{print}{MCS}(x, ...)
\method{show}{MCS}(x, ...)
\method{summary}{MCS}(object, ...)
\method{plot}{MCS}(x, xlab="p", ylab="MCS", ...)
\method{print}{bootMCS}(x, ...)
\method{show}{bootMCS}(x, ...)
\method{summary}{bootMCS}(object, alpha=0.05, ...)
\method{plot}{bootMCS}(x, xlab="p", ylab="MCS", alpha=0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ A matrix of numeric variables. }
  \item{p}{ The quantiles at which to evaluate. }
  \item{R}{The number of bootstrap samples to run. Defaults to
          \code{R = 100}.}
  \item{trace}{How often to inform the user of progress. Defaults to
          \code{trace = 10}.}
  \item{x, object}{An object of class \code{MCS} or \code{bootMCS}.}
  \item{xlab, ylab}{Axis labels.}
  \item{alpha}{A (1 - alpha)\% pointwise confidence interval will be produced.
               Defaults to \code{alpha = 0.05}.}
  \item{...}{Optional arguments to be passed into methods.}
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
plot(D)
# Following lines commented out to keep CRAN robots happy
#Dmcs <- bootMCS(D[, 5:6])
#Dmcs
#plot(Dmcs)
}
\keyword{ multivariate }

