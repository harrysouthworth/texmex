% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrl.R
\name{mrl}
\alias{mrl}
\alias{print.mrl}
\alias{summary.mrl}
\alias{print.summary.mrl}
\alias{plot.mrl}
\alias{ggplot.mrl}
\title{Mean residual life plot}
\usage{
mrl(data, umin = min(data), umax = max(data) - 0.1, nint = 100,
alpha=.050)
\method{print}{mrl}(x, ...)
\method{print}{summary.mrl}(x, ...)
\method{summary}{mrl}(object, ...)
\method{plot}{mrl}(x, xlab="Threshold", ylab="Mean excess", ...)
\method{ggplot}{mrl}(data, mapping, xlab = "Threshold",
  ylab = "Mean excess", main=NULL,fill="orange", col="blue",
  rug=TRUE, addNexcesses=TRUE, textsize=4, ..., environment)
}
\arguments{
\item{data}{A numeric vector.}

\item{umin}{The minimum value over which to threshold the data.}

\item{umax}{The maximum value over which to threshold the data.}

\item{nint}{The number of points at which to compute the plot.}

\item{alpha}{Used to determine coverage of confidence interval to plot.
Defaults to plotting a 95\% interval.}

\item{x, object}{Arguments to print, summary and plot functions.}

\item{xlab}{Label for the x-axis. Defaults to \code{xlab="Threshold"}.}

\item{ylab}{Label for the y-axis. Defaults to \code{ylab="Mean excess"}.}

\item{\dots}{Optional arguments to \code{plot}.}

\item{col}{Colour of the line on the MRL plot.}

\item{rug}{Whether to add raw data as a rug along axis of plot.}

\item{fill}{Colour of the pointwise confidence region on the MRL plot.}

\item{main}{Main title.}

\item{addNexcesses}{Whether to annotate the plot with the numbers of
excesses over increasing thresholds. Defaults to \code{addNexcesses=TRUE}.}

\item{textsize}{Size of text on the plot (ggplot). Defaults to
\code{textsize=4}.}

\item{mapping, environment}{Not used.}
}
\value{
A list with two components. \code{data} is the original data,
\code{mrl} is a matrix containing information to produce the mean residual
life plot.
}
\description{
Calculate mean residual life and plot it to aid the identification of a
threshold over which to fit a generalized Pareto distribution
}
\details{
Threshold choice for the fitting of the GPD is guided by the shape of the
Mean Residual Life plot.  A threshold which is suitably high will have a
corresponding mrl plot which is approximately linear in shape above the
threshold (up to sampling variation).
}
\note{
The function was originally written by Stuart Coles and appears in the
\code{ismev} package. This version modified by Harry Southworth to allow
more control over the appearance of the plot.
}
\references{
S. Coles, An Introduction to Statistical Modeling of Extreme
Values, Springer, 2001
}
\author{
Janet E. Heffernan, Harry Southworth
}
\keyword{models}
