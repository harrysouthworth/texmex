\name{plot.bgpd}
\alias{plot.bgpd}

\title{ Plot a \code{bgpd} object }
\description{
  Plots of posterior densities, and diagnostic plots, for an object of class
  \code{bgpd}.
}
\usage{
plot.bgpd(x, which.plots = 1:3, density.adjust = 2, print.seed = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class \code{bgpd}. }
  \item{which.plots}{ Specifiy which plots to create. Should be a numeric vector:
        1 indicates plots of posterior densities; 2 indicates plots of the traces
        of the Markov chains and their cumulative means; 3 indicates plots of the
        autocorrelation plots of the chains. }
  \item{density.adjust}{ A constant used to adjust the degree of smoothing done
        when computing the kernel density estimates. Defaults to 2. This argument
        is passed to \code{density}. See the help file for \code{density} for 
        more information. }
  \item{print.seed}{ Whether or not to print the seed that was used to the random
        number generator at the foot of each plot. Defaults to \code{FALSE}. }
  \item{\dots}{ For compliance with \code{plot}. Not used. }
}
\details{
  The posterior kernel density estimates are plotted, using more smoothing than
  \code{density} uses by default due to the Markov chain producing multiple
        identical values, resulting in a spikey kernel density estimate.
        
        If the chains have converged on the posterior distributions, the trace
        plots should look like `fat hairy caterpillars' and their cumulative 
        means should converge rapidly. Moreover, the autocorrelation functions
        should converge quickly to zero.
}
\value{ NULL }
\author{ Harry Southworth }
\seealso{ \code{\link{bgpd}}, \code{\link{summary.bgpd}} }
\examples{
x <- rnorm(1000)
mod <- bgpd(x, qu=0.7)
plot(mod)
}
\keyword{ models }

