\name{MCS}
\alias{MCS}
\alias{bootMCS}
\alias{plot.MCS}
\alias{summary.bootMCS}
\alias{plot.bootMCS}
\title{ Multivariate conditional Spearman's rho }
\description{
  Compute multivariate conditional Spearman's rho over a
  range of quantiles.
}
\usage{
MCS(X, p = seq(0.1, 0.9, by = 0.1))
bootMCS(X, p = seq(0.1, 0.9, by = 0.1), R = 100, trace = 10)

\method{plot}{MCS}(x, xlab="p", ylab="MCS", ...)

\method{summary}{bootMCS}(object, alpha=0.05, ...)
\method{plot}{bootMCS}(x, xlab="p", ylab="MCS", alpha=0.05, ylim, ...)
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
  \item{ylim}{Plotting limits for bootstrap plot.}
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
 MCS returns an object of class \code{MCS}.  There are plot and summary methods available for this class.
 
  \item{MCS }{The estimated correlations.}
  \item{p }{The quantiles at which the correlations were evaluated at}
  \item{call}{The function call used.}
  
  bootMCS returns an object of class \code{bootMCS}. There are plot and summary methods available for this class.
  
  \item{replicates}{Bootstrap replicates.}
  \item{p }{The quantiles at which the correlations were evaluated at}
  \item{R}{Number of bootstrap samples.}
  \item{call}{The function call used.}

   }
\references{ F. Schmid and R. Schmidt, Multivariate conditional versions
of Spearman's rho and related measures of tail dependence, Journal of
Multivariate Analysis, 98, 1123 -- 1140, 2007 }
\author{ Yiannis Papastathopoulos, Harry Southworth }

\seealso{ \code{\link{chi}} }
\examples{
D <- liver[liver$dose == "D",]
plot(D)
# Following lines commented out to keep CRAN happy
#Dmcs <- bootMCS(D[, 5:6])
#Dmcs
#plot(Dmcs)
}
\keyword{ multivariate }

