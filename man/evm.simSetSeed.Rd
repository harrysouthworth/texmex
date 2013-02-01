\name{evm.simSetSeed}
\alias{evm.simSetSeed}

\title{Set the seed from a fitted evm.sim object.}
\description{
Set the seed from a fitted \code{evm.sim} object to ensure reproducibility of output.
}
\usage{
evm.simSetSeed(x)
}

\arguments{
  \item{x}{An object of class \code{evm.sim}, as returned by \code{evm}
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
# evm.simSetSeed(mod)
# mod1 <- evm(data, qu=.7, method="simulate") # this produces the same MCMC output as mod
}
