\name{thinAndBurn}
\alias{thinAndBurn}
\alias{thinAndBurn.bgpd}
\title{
  Process Metropolis output from GPD fitting to discard unwanted observations.
}
\description{
  Process observations from Metropolis fitting of GPD models, to thin the output and discard observations from burn-in period.
}
\usage{
\method{thinAndBurn}{bgpd}(object, burn, thin)
}
\arguments{
  \item{object}{
Object of class 'bgpd' as returned by \code{gpd} called with \code{method="simulate"}.
}
  \item{thin}{
\code{thin} or its reciprocal must be a positive integer.  If integer valued, this specifies the frequency of observations from the simulated Markov Chain which will be retained.  If specified as a proportion, this is the proportion of values which will be retained. For no thinning use \code{thin=1}.
}
  \item{burn}{
The number of observations from the simulated Markov Chain to be discarded as burn-in. Must be a non-negative integer, for no burn-in use \code{burn=0}.
}
}

\value{
Object of class \code{bgpd}.  See Value returned by \code{\link{gpd}} using \code{method = "simulate"} for details.

Note that the original chain is not discarded when this function is called: \code{thinAndBurn} can be called recursively on the original object with different values of \code{burn} and \code{thin} without the this object getting progressively smaller!
}
\author{
Harry Southworth, Janet E. Heffernan 
}

\seealso{
\code{\link{gpd}}
}
\examples{
  x <- rnorm(1000)
  mod <- gpd(x, qu=.7, method="sim")
  mod
  par(mfrow=c(3, 2))
  plot(mod)
  mod1 <- thinAndBurn(mod,burn=1000, thin=5)
  plot(mod1)
}