\name{evmSimSetSeed}
\alias{evmSimSetSeed}

\title{Set the seed from a fitted evmSim object.}
\description{
Set the seed from a fitted \code{evmSim} object to ensure reproducibility of output.
}
\usage{
evmSimSetSeed(x)
}

\arguments{
  \item{x}{An object of class \code{evmSim}, as returned by \code{evm}
     using \code{method = "simulate"}.}
}
\details{
Sets the seed to the value used to fit the model.
}
\author{
Harry Southworth
}

\seealso{
\code{\link{evm}}
}
\examples{
# data <- rnorm(1000)
# mod <- evm(data, qu=.7, method="simulate")
# evmSimSetSeed(mod)
# mod1 <- evm(data, qu=.7, method="simulate") # this produces the same MCMC output as mod
}
