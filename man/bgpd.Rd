\name{bgpd}
\alias{bgpd}
\alias{print.bgpd}
\alias{summary.bgpd}
\alias{plot.bgpd}
\alias{coef.bgpd}
\alias{coefficients.bgpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayesian inference for the generalized Pareto distribution }
\description{
  Simulate posterior distributions for the parameters in a generalized
  Pareto distribution, optionally with explanatory variables
}
\usage{
bgpd(y, data, th, qu, phi = ~1, xi = ~1, prior = "gaussian", priorParameters = NULL,
     iter = 10500, burn = 500, thin = 1, jump.const, start, trace = 1000, verbose=TRUE)
\method{print}{bgpd}(x, print.seed=FALSE, ...)
\method{summary}{bgpd}(object, ...)
\method{plot}{bgpd}(x, which.plots=1:3, density.adjust=2, print.seed=FALSE, ...)
\method{coef}{bgpd}(object, ...)
\method{coefficients}{bgpd}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ The data vector, or the name of a variable in \code{data}. }
  \item{data}{ A data object containing the data to be modelled and any covariates. }
  \item{th}{ The threshold to be applied to \code{y}. }
  \item{qu}{An alternative to \code{th}, specifying the quantile of the \code{y}
            above which the model should be fit.}
  \item{phi}{ A formula for modelling the log of the scale parameter in the generalized Pareto distribution.
              Working with the log was found to improve stability of computations
              when the sample size is small, and makes Gaussian priors sensible.
              Defaults to \code{phi = ~ 1} - i.e. no covariates.}
  \item{xi}{ A formula for modelling the log of the scale parameter in the generalized Pareto distribution.
                Defaults to \code{phi = ~ 1} - i.e. no covariates. }
  \item{prior}{ The type of prior distribution to use. At present, only ``none''
                 or ``gaussian'' (the default) is accepted. }
  \item{priorParameters}{ A list with two components. The first should be 
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
  \item{trace}{ The size of the interval of numbers of iterations at which to
                inform the user of progress. Defaults to \code{trace=1000}. }
  \item{verbose}{Whether or not to print progress to screen. Defaults to
                \code{verbose=TRUE}.}
  \item{print.seed}{Whether or not to print the seed used in the simulations,
                    or to annotate the plots with it. Defaults to
                    \code{print.seed=FALSE}.}
  \item{which.plots}{Which plots to produce. Option 1 gives kernel density estimates,
                     2 gives traces of the Markov chains with superimposed
                     cumulative means, 3 gives autocorrelation functions.
                     Defaults to \code{which.plots=1:3}.}
  \item{density.adjust}{Passed into \code{density}. Controls the amount of
                     smoothing of the kernel density estimate. Defaults to
                     \code{density.adjust=2}.}

}
\details{
  If \code{start} is not provided, the maximum penalized likelhood point
  estimates are computed and used. The simulation is done by a Metropolis
  algorithm.

  When plotting the object, f the chains have converged on the posterior distributions, the trace
        plots should look like `fat hairy caterpillars' and their cumulative 
        means should converge rapidly. Moreover, the autocorrelation functions
        should converge quickly to zero.

  When printing or summarizing the object, 
  posterior means and standard deviations are computed. Posterior means
are also returned by the \code{coef} method. Depending on what you
  want to do and what the posterior distributions look like (try using \code{plot.bgpd})
  you might want to work with quantiles of the poseterior distributions instead.

}
\value{
    \item{call}{The call to \code{bgpd} that produced the object.}
    \item{threshold}{The threshold above which the model was fit.}
    \item{map}{The point estimates found by maximum penalized likelihood
               and which were used as the starting point for the Markov
               chain.}
    \item{burn}{The number of steps of the Markov chain that are to be 
                treated as the burn-in and not used in inferences.}
    \item{thin}{The degree of thinning used.}
    \item{chains}{The entire Markov chain generated by the Metropolis
                  algorithm.}
    \item{param}{The remainder of the chain after deleting the burn-in
                 and applying any thinning.}
    \item{X.phi}{The design matrix for the log of the scale parameter.}
    \item{X.xi}{The design matrix for the log of the scale parameter.}
    \item{acceptance}{The proportion of proposals that were accepted by
                      the Metropolis algorithm.}
    \item{seed}{The seed used by the random number generator.}
}
\author{ Harry Southworth }
\seealso{ \code{\link{gpd}} }
\examples{
  x <- rnorm(1000)
  mod <- bgpd(x, qu=0.7)
  plot(mod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

