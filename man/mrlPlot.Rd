\name{mrlPlot}
\alias{mrlPlot}
\title{ Mean residual life plot }
\description{
  Produce a mean residual life plot to aid the identification of a threshold over which
  to fit a generalized Pareto distribution
}
\usage{
mrlPlot(data, umin = min(data), umax = max(data) - 0.1, conf = 0.95, nint = 100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A numeric vector. }
  \item{umin}{ The minimum value over which to threshold the data. }
  \item{umax}{ The maximum value over which to threshold the data. }
  \item{conf}{ The coverage of confidence interval to plot. Defaults to
         plotting a 95\% interval. }
  \item{nint}{ The number of points at which to compute the plot. }
  \item{\dots}{ Optional arguments to \code{plot} }
}
\details{
  Threshold choice for the fitting of the GPD is guided by the shape of the Mean Residual Life plot.  A threshold which is suitably high will have a corresponding mrl plot which is approximately linear in shape above the threshold (up to sampling variation).

  The function was originally written by Janet E. Heffernan and appears in
  the \code{ismev} package. This version modified by Harry Southworth to
  allow more control over the appearance of the plot.
}
\value{
NULL
}
\author{ Janet E. Heffernan, Harry Southworth }

\keyword{ models }

