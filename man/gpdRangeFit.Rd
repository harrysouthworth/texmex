\name{gpdRangeFit}
\alias{gpdRangeFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate generalized Pareto distribution parameters over a range of values }
\description{
  Estimate generalized Pareto distribution parameters over a range of values, using
  maximum (penalized) likelihood.
}
\usage{
gpdRangeFit(data, umin, umax, nint = 10, show = FALSE, xlab = "Threshold", ylab = NULL, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The data vector to be modelled. }
  \item{umin}{ The minimum threshold to esimate the parameters over. }
  \item{umax}{ The maximum threshold to esimate the parameters over. }
  \item{nint}{ The number of points at which to perform the estimation. }
  \item{show}{ ~~Describe \code{show} here~~ }
  \item{xlab}{ Label for the x-axis. }
  \item{ylab}{ Label for the y-axis. }
  \item{main}{ The main title. }
  \item{\dots}{ Arguments to \code{plot} }
}
\details{
  This is Janet E Heffernan's \code{gpd.fitrange}, as it appears in the
  \code{ismev} package, with some fairly minor changes. The function uses
  \code{gpd} internally and uses the default options for that function.
  
}

\author{ Janet E Heffernan, Harry Southworth }
\seealso{ gpd }
\examples{

}
\keyword{ models }

