% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot.boot.R
\name{ggplot.evmBoot}
\alias{ggplot.evmBoot}
\alias{ggbootdensplots}
\title{Diagnostic plots for the replicate estimated parameter values in an evmBoot object}
\usage{
\method{ggplot}{evmBoot}(
  data = NULL,
  mapping,
  denscol = "light blue",
  histcol = "dark blue",
  linecol = "orange",
  plot.it = TRUE,
  ...,
  environment
)
}
\arguments{
\item{data}{An object of class 'evmBoot'.}

\item{mapping, environment}{ignored}

\item{denscol}{Colour for the densities. Defaults to 'light blue'.}

\item{histcol}{Colour for the histograms. Defaults to 'dark blue'.}

\item{linecol}{Colour for the point estimate lines. Decaults to 'orange'.}

\item{plot.it}{Whether or not to actually print the plots. Defaults
to \code{plot.it=TRUE}.  If \code{plot.it=FALSE}, you might
want to control the layout. Do this with
\code{do.call("grid.arrange", c(plots, ncol=2))}, for example,
where \code{plots} is the objected returned by
\code{ggplot.evmBoot}.}

\item{...}{Additional arguments to \code{ggplot}, currently unused.}
}
\description{
Diagnostic plots for the replicate estimated parameter values in an evmBoot object
}
\keyword{hplot}
