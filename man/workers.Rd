\name{texmexWorkers}
\alias{texmexMethod}
\alias{texmexPrior}
\alias{texmexTrace}
\alias{texmexPrepareData}
\alias{texmexThresholdData}
\alias{texmexPriorParameters}
\alias{findFormulae}
\alias{texmexParameters}
\alias{texmexGetParam}
\title{Worker functions for texmex}
\description{Worker functions not intended for direct use.}
\usage{
texmexMethod(method)
texmexPrior(prior, penalty, method, pp)
texmexTrace(trace, method)
texmexPrepareData(y, data, params)
texmexThresholdData(threshold, data)
texmexPriorParameters(prior, priorParameters, data)
findFormulae(call, ...)
texmexParameters(call, fam, ...)
texmexGetParam(data, co)
}

\arguments{
\item{method}{Character string describing the method to use. The
  function attempts to interpret it and returns "o", "s" or "b".}
\item{prior}{Character string describing the type of prior distribution
  to use. The function attempts to interpret it and returns a
  standardized version.}
\item{penalty}{Character string describing the type of penalty to use.}
\item{pp}{The \code{priorParameters} object.}
\item{trace}{Information on whether to report progress to the user.}
\item{y}{Vector of data to be modelled.}
\item{data}{A \code{data.frame} containing covariates. Otherwise a list
  of \code{data.frame}s already preprocessed.}
\item{params}{Names of the parameters in the model.}
\item{threshold}{The threshold above which to model the data.}
\item{priorParameters}{A \code{list} containing information on the prior
  distribution.}
\item{call}{An object of class "call".}
\item{fam}{An object of class "texmexFamily".}
\item{co}{Model coefficients.}
\item{...}{Other arguments being passed between functions.}
}
\details{These are functions used for processing user input and data
  prior to calling \code{evm.fit} and are not intended for end-user usage.
}
\value{Processed data and other information. }
\author{ Harry Southworth }

\keyword{multivariate}

