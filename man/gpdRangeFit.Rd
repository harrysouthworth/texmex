\name{gpdRangeFit}
\alias{gpdRangeFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate generalized Pareto distribution parameters over a range of values }
\description{
  Estimate generalized Pareto distribution parameters over a range of values, using
  maximum (penalized) likelihood.
}
\usage{
gpdRangeFit(data, umin, umax, nint = 10, 
            penalty = "gaussian", priorParameters = NULL, alpha=0.05,
xlab = "Threshold", ylab = NULL, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The data vector to be modelled. }
  \item{umin}{ The minimum threshold above which to esimate the parameters. }
  \item{umax}{ The maximum threshold above which to esimate the parameters. }
  \item{nint}{ The number of thresholds at which to perform the estimation. }
  \item{penalty}{The type of penalty to be used in the maximum penalized
               likelihood estimation. Should be either ``gaussian'' or 
               ``none''. Defaults to ``gaussian''.}
  \item{priorParameters}{Parameters to be used for the penalty function.
                 See the help for \code{\link{gpd}} for more informaiton.}
  \item{alpha}{(1 - alpha)\% confidence intervals will be plotted with the 
               point estimates. Defaults to \code{alpha = 0.05}.}
  \item{xlab}{ Label for the x-axis. }
  \item{ylab}{ Label for the y-axis. }
  \item{main}{ The main title. }
  \item{\dots}{ Arguments to \code{plot} }
}
\details{
  This is Stuart Coles's \code{gpd.fitrange}, as it appears in the
  \code{ismev} package, with some fairly minor changes. The function uses
  \code{gpd} internally and uses the default options for that function.
  
}

\author{ Stuart Coles, Janet E Heffernan, Harry Southworth }
\seealso{ gpd }
\examples{
par(mfrow=c(1,2))
gpdRangeFit(rain)
}
\keyword{ models }

