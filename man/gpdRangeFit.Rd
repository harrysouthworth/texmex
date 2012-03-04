\name{gpdRangeFit}
\alias{gpdRangeFit}
\title{ Estimate generalized Pareto distribution parameters over a range of values }
\description{
  Estimate generalized Pareto distribution parameters over a range of values, using
  maximum (penalized) likelihood.
}
\usage{
gpdRangeFit(data, umin=quantile(data, .05), umax=quantile(data, .95), nint = 10, 
            penalty = "gaussian", priorParameters = NULL, alpha=0.05,
            xlab = "Threshold", ylab = NULL, main = NULL, addNexcesses=TRUE, ...)
}

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
  \item{addNexcesses}{ Annotate top axis with numbers of threshold excesses arising with the corresponding values of threshold on the bottom axis.}
  \item{\dots}{ Arguments to \code{plot} }
}
\details{
  This is Stuart Coles' \code{gpd.fitrange}, as it appears in the
  \code{ismev} package, with some fairly minor changes. The function uses
  \code{gpd} internally and uses the default options for that function.
  
  Note this function does not extend to assessing model fit when there are covariates included in the mdoel.
}

\author{ Stuart Coles, Janet E Heffernan, Harry Southworth }
\seealso{ \code{\link{gpd}} }
\examples{
par(mfrow=c(1,2))
gpdRangeFit(rain)
}
\keyword{ models }

