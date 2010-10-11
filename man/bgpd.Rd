\name{bgpd}
\alias{bgpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayesian inference for the generalized Pareto distribution }
\description{
  Simulate posterior distributions for the parameters in a generalized
  Pareto distribution, optionally with explanatory variables
}
\usage{
bgpd(y, data, th, phi = ~1, xi = ~1, prior = "gaussian", gaussianParameters = NULL, iter = 10500, burn = 500, thin = 1, jump.const, start, trace = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ The data vector, or the name of a variable in \code{data}. }
  \item{data}{ A data object containing the data to be modelled and any covariates. }
  \item{th}{ The threshold to be applied to \code{y}. }
  \item{phi}{ A formula for modelling the log of the scale parameter in the generalized Pareto distribution.
              Working with the log was found to improve stability of computations
              when the sample size is small, and makes Gaussian priors sensible.
              Defaults to \code{phi = ~ 1} - i.e. no covariates.}
  \item{xi}{ A formula for modelling the log of the scale parameter in the generalized Pareto distribution.
                Defaults to \code{phi = ~ 1} - i.e. no covariates. }
  \item{prior}{ The type of prior distribution to use. At present, only ``none''
                 or ``gaussian'' (the default) is accepted. }
  \item{gaussianParameters}{ A list with two components. The first should be 
                a vector of means, the second should be a covariance matrix.
                These represent the parameters in the Gaussian prior distribution.
                If not supplied: all default prior means are zero;
                all default prior variances are $10^4$;all covariances
                are zero. }
  \item{iter}{ The number of steps for the Metropolis algorith to take. }
  \item{burn}{ The number of initial steps to be discarded. }
  \item{thin}{ The degree of thinning of the resulting Markov chains. Defaults
               to 1 (no thinning). A value of 0.5 (for example) would result in
               every other step being discarded. }
  \item{jump.const}{ Control parameter for the Metropolis algorithm. }
  \item{start}{ Starting values for the parameters. }
  \item{trace}{ Whether or not to print progress to screen. }
}
\details{
  If \code{start} is not provided, the maximum penalized likelhood point
  estimates are computed and used. The simulation is done by a Metropolis
  algorithm.
}
\value{
  A list with class ``bgbm''. There are print, summary and plot functions
  available.
}
\author{ Harry Southworth }
\seealso{ \code{\link{gpd}} }
\examples{
  x <- rnorm(1000)
  mod <- bgpd(x, th=quantile(x, .7))
  plot(mod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

