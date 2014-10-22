\name{simulate.evmOpt}
\alias{simulate.evmOpt}
\alias{simulate.evmSim}
\alias{simulate.evmBoot}
\title{Simulate from a fitted evm object}
\description{
  Simulate random numbers from a fitted evm object
}
\usage{
\method{simulate}{evmOpt}(object, nsim=1, seed=NULL, param=NULL, ...)
\method{simulate}{evmSim}(object, nsim=NULL, seed=NULL, ...)
\method{simulate}{evmBoot}(object, nsim=NULL, seed=NULL, ...)
}

\arguments{
  \item{object}{A fitted evm object having class 'evmOpt', 'evmSim' or 'evmBoot'.}
  \item{nsim}{The number of random numbers to simulate. For \code{simulate.evmOpt},
        defaults to \code{nsim=1}. For \code{simulate.evmSim} and \code{simulate.evmBoot},
        it is ignored and the number of simulated values is inferred from the number of
        simulated parameters left in the Markov chains after removing the burn-in and
        thinning, or the number of bootstrap replicates, respectively.}
  \item{seed}{An integer to be passed to \code{set.seed}. Defaults to \code{seed=NULL}.}
  \item{param}{Parameters to use in the random number generator. Defaults to
        \code{param=NULL} in which case the parameters from the fitted model are used.
        For \code{simulate.evmSim} and \code{simulate.evmBoot}, this argument is not
        available and the simulated parameters or replicates are used.}
  \item{...}{Unused.}
}
\value{
  A vector or random numbers simulated from the fitted model object.
}

\author{Paul Metcalfe, Harry Southworth}
\seealso{ \code{\link{evm}} }
\examples{
mod <- evm(rain, qu=.95)
hist(simulate(mod, 100))
}
\keyword{models}

