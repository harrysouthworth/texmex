\name{methods}
\alias{coef.evmSim}
\alias{coef.evmBoot}
\alias{coefficients.evmBoot}
\alias{coef.mex}
\alias{coef.migpd}
\alias{coefficients.evmSim}
\alias{coefficients.evmOpt}
\alias{coefficients.mex}
\alias{coefficients.migpd}

\alias{plot.extremalIndex}
\alias{plot.ppevm}
\alias{plot.qqevm}
\alias{plot.mexRangeFit}
\alias{plot.hist.evmOpt}

\alias{print.MCS}
\alias{print.bootMCS}
\alias{print.chi}
\alias{print.copula}
\alias{print.declustered}
\alias{print.extremalIndex}
\alias{print.mex}
\alias{print.mexDependence}
\alias{print.migpd}
\alias{print.predict.mex}
\alias{print.summary.evmSim}
\alias{print.summary.evmOpt}
\alias{print.summary.evmBoot}
\alias{print.summary.predict.mex}
\alias{print.hist.evmOpt}
\alias{print.ppevm}

\alias{summary.MCS}
\alias{summary.copula}
\alias{summary.mex}
\alias{summary.migpd}

\title{ Methods for texmex objects}
\description{
  Methods for texmex objects.
}
\usage{

\method{coef}{evmSim}(object,...)
\method{coef}{mex}(object,...)
\method{coef}{migpd}(object,...)
\method{coefficients}{evmSim}(object,...)
\method{coefficients}{evmOpt}(object,...)
\method{coefficients}{mex}(object,...)
\method{coefficients}{migpd}(object,...)

\method{plot}{extremalIndex}(x,...)
\method{plot}{mexRangeFit}(x, col=2, bootcol="grey", addNexcesses=TRUE, ...)

\method{print}{MCS}(x,...)
\method{print}{bootMCS}(x,...)
\method{print}{chi}(x,...)
\method{print}{copula}(x,...)
\method{print}{declustered}(x,...)
\method{print}{extremalIndex}(x,...)
\method{print}{mex}(x,...)
\method{print}{mexDependence}(x,...)
\method{print}{migpd}(x,...)
\method{print}{predict.mex}(x,...)
\method{print}{summary.evmSim}(x,...)
\method{print}{summary.evmOpt}(x,digits=3,...)
\method{print}{summary.predict.mex}(x,...)
\method{summary}{MCS}(object,...)
\method{summary}{copula}(object,...)
\method{summary}{mex}(object,...)
\method{summary}{migpd}(object,verbose=TRUE,...)
}

\arguments{
  \item{x}{ A texmex object of class indicated by method type. }
  \item{object}{ A texmex object of class indicated by method type.}
  \item{digits}{ Number of digits for printing.}
  \item{verbose}{ Logical, whether to show model fitting summary in addition to just model parameter estimates.}
  \item{col}{ The colour of points on the plot (the point estimates in \code{plot.mexRangeFit}).}
  \item{bootcol}{ The colour of points for point estimates from bootstrap samples in \code{plot.mexRangeFit}.}
  \item{addNexcesses}{Whether or not to print the number of threshold excesses in \code{plot.mexRangeFit}.
        Defaults to \code{TRUE}.}
  \item{...}{Further arguments to methods.}
  }

\author{ Harry Southworth, Janet E. Heffernan }