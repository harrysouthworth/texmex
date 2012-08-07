\name{validate.texmex}
\alias{validate.texmex}
\alias{test.gpd}
\alias{test.mex}
\alias{test.bootgpd}
\alias{test.bootmex}
\alias{test.endPoint}
\alias{test.extremalIndex}
\alias{test.gpdRangeFit}
\alias{test.MCS}
\alias{test.mexDependence}
\alias{test.mexRangeFit}
\alias{test.migpd}
\alias{test.migpdCoefs}
\alias{test.mrlPlot}
\alias{test.pgpd}
\alias{test.plot.bootmex}
\alias{test.plot.chi}
\alias{test.plot.gpd}
\alias{test.plot.mexDependence}
\alias{test.plot.migpd}
\alias{test.plot.lp.gpd}
\alias{test.plot.predict.mex}
\alias{test.plotrl.gpd}
\alias{test.predict.bgpd}
\alias{test.predict.bootgpd}
\alias{test.predict.gpd}
\alias{test.predict.mex}
\alias{test.qgpd}
\alias{test.rgpd}
\alias{test.thinAndBurn.bgpd}
\alias{test.chi}
\alias{test.copula}
\alias{test.dgpd}
\alias{test.info.gpd}
\alias{test.plot.mex}
\alias{test.revTransform}
\alias{test.predict.gpd}
\alias{test.predict.bgpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Validate the texmex package
}
\description{
Using package \code{RUnit}, run tests to compare behaviour of \code{texmex} with
expected behaviour
}
\usage{
validate.texmex()
}
\details{
The function uses functionality in the \code{RUnit} package to perform many
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
original published versions to further check that the \code{texmex} reproduces 
expected output.

If all is well, the result should indicate that no errors or fails have
occurred.
}
\value{
An object of class \code{RUnitTestData}. See the help for \code{RUnit} for
details.

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
Around 400 checks are performed, some of which run MCMC or bootstrap algorithms.
As a result, the function takes a while to run.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{texmex-package}}
}
\examples{
# Uncomment the following lines to run - commented out to make CRAN happy
## open device to save plotted output
#library(RUnit)
#pdf("texmexPlots.pdf")
#val <- validate.texmex()
#printHTMLProtocol(val, "texmexReport.html")
## finish writing to device
#dev.off()
}
\keyword{ models }
\keyword{multivariate}

