\name{mrl}
\alias{mrl}
\alias{print.mrl}
\alias{summary.mrl}
\alias{plot.mrl}
\title{ Mean residual life plot }
\description{
  Calculate mean residual life and plot it to aid the
  identification of a threshold over which to fit a
  generalized Pareto distribution
}
\usage{
mrl(data, umin = min(data), umax = max(data) - 0.1, nint = 100, alpha=.050)
\method{print}{mrl}(x, ...)
\method{summary}{mrl}(object, ...)
\method{plot}{mrl}(x, xlab="Threshold", ylab="Mean excess", ...)
}

\arguments{
  \item{data}{ A numeric vector. }
  \item{umin}{ The minimum value over which to threshold the data. }
  \item{umax}{ The maximum value over which to threshold the data. }
  \item{alpha}{ Used to determine coverage of confidence interval to plot. Defaults to
         plotting a 95\% interval. }
  \item{nint}{ The number of points at which to compute the plot. }
  \item{x, object}{Arguments to print, summary and plot functions.}
  \item{xlab}{Label for the x-axis. Defaults to \code{xlab="Threshold"}.}
  \item{ylab}{Label for the y-axis. Defaults to \code{ylab="Mean excess"}.}
  \item{\dots}{ Optional arguments to \code{plot} }
}
\details{
  Threshold choice for the fitting of the GPD is guided by the shape of the Mean Residual Life plot.  A threshold which is suitably high will have a corresponding mrl plot which is approximately linear in shape above the threshold (up to sampling variation).
}
\value{
A list with two components. \code{data} is the original data,
\code{mrl} is a matrix containing information to produce the
mean residual life plot.
}
\note{  The function was originally written by Stuart Coles and appears in
  the \code{ismev} package. This version modified by Harry Southworth to
  allow more control over the appearance of the plot.}
\author{ Janet E. Heffernan, Harry Southworth }
\references{
S. Coles, An Introduction to Statistical Modeling of Extreme Values, Springer, 2001
}
\keyword{ models }

