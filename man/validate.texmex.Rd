\name{validate.texmex}
\alias{validate.texmex}
\alias{test.evm}
\alias{test.mex}
\alias{test.evmBoot}
\alias{test.bootmex}
\alias{test.endPoint}
\alias{test.extremalIndex}
\alias{test.gpdRangeFit}
\alias{test.MCS}
\alias{test.mexDependence}
\alias{test.mexRangeFit}
\alias{test.migpd}
\alias{test.migpdCoefs}
\alias{test.mrl}
\alias{test.pgpd}
\alias{test.plot.bootmex}
\alias{test.plot.chi}
\alias{test.plot.gpd}
\alias{test.plot.mexDependence}
\alias{test.plot.migpd}
\alias{test.plot.lp.gpd}
\alias{test.plot.predict.mex}
\alias{test.plotrl.evmOpt}

\alias{test.evmOpt}
\alias{test.plot.evmOpt}
\alias{test.plot.lp.evmOpt}
\alias{test.plotrl.evm}
\alias{test.predict.evm}
\alias{test.thinAndBurn.evmSim}

\alias{test.predict.evmSim}
\alias{test.predict.evmBoot}
\alias{test.predict.evmOpt}
\alias{test.predict.mex}
\alias{test.qgpd}
\alias{test.rgpd}
\alias{test.thinAndBurn.bgpd}
\alias{test.chi}
\alias{test.copula}
\alias{test.dgpd}
\alias{test.gpd.info}
\alias{test.plot.mex}
\alias{test.revTransform}
\alias{test.predict.evmOpt}
\alias{test.predict.evmSim}

\alias{test.closures}
\alias{test.log1prel}
\alias{test.pgev}
\alias{test.qgev}
\alias{test.rgev}
\alias{test.specfun.safe.product}

\alias{test.exprel}
\alias{test.log1prel}
\alias{test.specfun.safe.product}

\alias{test.plot.predict.evm}

\title{Validate the texmex package}
\description{
Using package \code{RUnit}, run tests to compare behaviour of \code{texmex} with
expected behaviour
}
\usage{
validate.texmex()
}
\details{
The function uses functionality in the \code{testthat} package to perform many
tests to ensure that \code{texmex} behaves as expected. These include tests
to make sure that \code{texmex} can reproduce published results (often parameter
estimates and model summaries that appear in Coles (2001) or in Heffernan
and Tawn (2004)), logical tests (for
example, to ensure that penalizing the likelihood has the expected effect on
parameter estimates), tests of the structure of objects, tests of proper
vectorization, tests to ensure that simulated posterior means are close to
MAP estimates, and tests to ensure that the code gives the same results as
independently coded calculations (sometimes using code that appears in the
\code{ismev} and \code{evd} packages, sometimes using code written by a third party, without reference to the \code{texmex} version).  A further side 
effect is the production of plots, many of which can be compared with their 
original published versions to further check that \code{texmex} reproduces 
expected output.

If all is well, the result should indicate that no errors or fails have
occurred.
}
\value{
A data frame of the summary test results.

Plots are also produced on the active plotting device; those plots 
that reproduce figures in the referenced material are labelled with the Figure 
numbers of the plot which they emulate.
}
\references{
S. Coles, An Introduction to the Statistical Modeling of Extreme Values, Springer, 2001

J. E. Heffernan and J. A. Tawn, A conditional approach for multivariate extreme
values (with discussion). Journal of the Royal Statistical Society, Series B, volume 66, 2004
}
\author{
Janet E. Heffernan, Harry Southworth
}
\note{
Around 1000 checks are performed, some of which run MCMC or bootstrap algorithms.
As a result, the function takes a while to run.
}

\seealso{
\code{\link{texmex-package}}
}
\examples{
# Uncomment the following lines to run
## open device to save plotted output
#pdf("texmexPlots.pdf")
#val <- validate.texmex()
#val
## finish writing to device
#dev.off()
}
\keyword{models}
\keyword{multivariate}

