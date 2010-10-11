\name{plot.gpd}
\alias{plot.gpd}
\alias{hist.gpd}
\alias{rl.gpd}
\alias{ppgpd}
\alias{qqgpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots for a `gpd' object }
\description{
  Diagnostic plots for a fitted generalized Pareto distribution
}
\usage{
plot.gpd(object, main = rep(NULL, 4), xlab = rep(NULL, 4))
}

\arguments{
  \item{object}{ An object of class `gpd' as returned by \code{gpd} }
  \item{main}{ Title for the plot }
  \item{xlab}{ X-axis label }
}
\details{
  Diagnostic plots for generalized Pareto models. A PP-plot, QQ-plot,
  histogram with superimposed generalized Pareto density estimate,
  and a return level plot with confidence interval are produced.
  
  The PP-plot and QQ-plot contain simulated pointwise confidence regions.
  By default, the region is a 95\% region based on 1000 simulated
  samples.
}
\value{
  NULL, but \code{qqplot} can be used to return a summary of
  the simulated confidence region.
}
\references{ S. Coles, An Introduction to the Statistical Modelling
     of Extreme Values, Springer, 2001 }
\author{ Harry Southworth, Janet E Heffernan }
\note{ Most of the code was originally written by Janet E Heffernan and
       appears in the \code{ismev} package.
}
\seealso{ \code{\link{gpd}} }
\examples{
  alt <- log(liver[liver$dose == "D", "ALT_M"])
  mod <- gpd(alt, th=quantile(alt, .7))
  par(mfrow=c(2, 2))
  plot(mod)
  args(qqgpd)
  myenvelope <- qqgpd(mod, plot=FALSE, nsim=10000, alpha=.1)
  names(myenvelope)
}
\keyword{ models }

