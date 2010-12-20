\name{texmex-internal}
\alias{hist.gpd}
\alias{rl.gpd}
\alias{qqgpd}
\alias{ppgpd}
\alias{u2gpd}
\alias{revGumbel}
\alias{mexGumbel}
\alias{gpd.fit}
\title{Internal functions for \code{texmex}}
\description{Internal functions used by the \code{texmex} package.}

\usage{
hist.gpd(x, xlab, ylab, main, ...)
rl.gpd(object, alpha = 0.05, xlab, ylab, main, pch=1, col=2, cex=0.75,
       linecol=4, cicol=0, polycol=15, smooth=TRUE)
qqgpd(object, nsim = 1000, alpha = 0.05, xlab, ylab, main, plot = TRUE,
     ylim = "auto", pch=1, col = 2, cex = 0.75, linecol = 4, cicol = 0,
     polycol = 15, smooth = TRUE)
ppgpd(object, nsim = 1000, alpha = 0.05, xlab, ylab, main,
     pch=1, col = 2, cex = 0.75, linecol = 4, cicol = 0,
     polycol = 15, smooth = TRUE)
u2gpd(u, p=1, th=0, sigma, xi)
mexGumbel(x, method = "mixture", divisor = "n+1", na.rm=TRUE)
revGumbel(x, data, qu, th=0, sigma=1, xi=0, method="mixture")
gpd.fit(y, th, X.phi, X.xi, penalty="none", start=NULL,
        priorParameters = NULL, maxit = 10000, trace = 0)
}

\arguments{
\item{x, object, data}{Object to be used by plot functions, vector to be converted.}
\item{xlab, ylab, ylim, main, pch, col, cex, linecol, cicol, polycol, ...}{Arguments to plot functions.}
\item{smooth}{Whether or not to smooth simulated confidence intervlas. 
    Defaults to \code{smooth = TRUE}.}
\item{alpha}{Control nominal coverage of condfidence intervals. Defaults to \code{alpha = 0.05}.}
\item{nsim}{Number of simulated datasets to use in computing confidence intervals.}
\item{plot}{Whether or not to plot the result.}
\item{u}{Uniform numbers to be converted to GPD numbers.}
\item{p, th, qu, sigma, xi}{Paramters of GPD distribution.}
\item{method}{How to convert. When \code{method = 'mixture'}, the upper tail of the
    distribution is modelled using a generalized Pareto distribution and the remainder
    is approximated using the empirical distribution. }
\item{divisor}{Divisor used in estimation of empirical distribution.}
\item{na.rm}{Whether or not to remove missing values.}
\item{y, X.phi, X.xi, penalty, start, priorParameters, maxit, trace}{Arguments
    supplied to \code{gpd}, \code{migpd} or \code{mex}, or inferred from
    those functions after some preprocessing.}
}
\details{
None of these functions are intended to be used explicitly.

The plotting functions are used internally by \code{\link{plot.gpd}}.

Some of the code is based on code that appears in the \code{ismev} package,
originally written by Stuart Coles.
}
\author{Harry Southworth, Janet E. Heffernan}



