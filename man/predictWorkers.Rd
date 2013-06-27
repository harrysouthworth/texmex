\name{predictWorkers}
\alias{addCov}
\alias{texmexMakeParams}
\alias{texmexMakeCovariance}
\alias{texmexPredictSE}
\alias{texmexMakeCI}
\alias{texmexMakeNewdataD}
\alias{texmexMakeCISim}
\title{Worker functions for prediction}
\description{Worker functions not intended for direct use.}
\usage{
addCov(res, X)
texmexMakeParams(co, data)
texmexMakeCovariance(cov, data)
texmexPredictSE(x)
texmexMakeCI(params, ses, alpha)
texmexMakeNewdataD(x, newdata)
texmexMakeCISim(x, alpha, object, sumfun, M)
}

\arguments{
  \item{res}{The actual predictions to have covariates added.}
  \item{X}{A design matrix.}
  \item{co}{Model coefficients.}
  \item{data}{A list of design matrices.}
  \item{cov}{Full covariance matrix for all parameters.}
  \item{x}{A covariance matrix. Otherwise an "evm" object.}
  \item{params}{Matrix of parameter values.}
  \item{ses}{Standard errors.}
  \item{alpha}{100(1-\code{alpha})\% confidence intervals will be produced.}
  \item{newdata}{Data for which predictions are required.}
  \item{object}{An "evm" object.}
  \item{sumfun}{A function for summarizing data.}
  \item{M}{The return level required.}
}
\details{These are functions used for processing information by
  \code{predict.evm} and are not intended for end-user usage.
}
\value{Processed data and other information.}
\author{Harry Southworth}

\keyword{models}

