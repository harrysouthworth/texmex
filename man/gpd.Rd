\name{gpd}
\alias{gpd}
\alias{print.gpd}
\alias{show.gpd}
\alias{summary.gpd}
\alias{plot.gpd}
\alias{coef.gpd}
\alias{AIC.gpd}
\alias{coefficients.gpd}
\alias{print.bgpd}
\alias{summary.bgpd}
\alias{plot.bgpd}
\alias{coef.bgpd}
\alias{coefficients.bgpd}
\alias{print.summary.bgpd}
\alias{print.summary.gpd}
\alias{show.summary.gpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Pareto distribution modelling }
\description{Likelihood based modelling and inference for the generalized
Pareto distribution, possibly with explanatory variables.
}
\usage{
gpd(y, data, th, qu, phi = ~1, xi = ~1, penalty = "gaussian",
    prior = "gaussian", method = "optimize", start = NULL, priorParameters = NULL,
    maxit = 10000, trace = NULL,
    iter = 10500, burn = 500, thin = 1, jump.const, verbose = TRUE)
\method{print}{gpd}(x, digits=max(3, getOption("digits") - 3), ...)
\method{show}{gpd}(x, digits=max(3, getOption("digits") - 3), ...)
\method{summary}{gpd}(object, nsim=1000, alpha=0.05, ...)
\method{coef}{gpd}(object, ...)
\method{plot}{gpd}(x, main=rep(NULL, 4), xlab=rep(NULL, 4), nsim=1000, alpha=0.05, ...)
\method{AIC}{gpd}(object, ..., k=2)
\method{coefficients}{gpd}(object, ...)
\method{print}{bgpd}(x, print.seed=FALSE, ...)
\method{summary}{bgpd}(object, ...)
\method{plot}{bgpd}(x, which.plots=1:3, density.adjust=2, print.seed=FALSE, ...)
\method{coef}{bgpd}(object, ...)
\method{coefficients}{bgpd}(object, ...)
\method{print}{summary.gpd}(x, dig=3, ...)
\method{print}{summary.bgpd}(x, ...)
\method{show}{summary.gpd}(x, dig=3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Either a numeric vector or the name of a variable in \code{data}. }
  \item{data}{ A data object containing \code{y} and any covariates. }
  \item{th}{ The threshold for \code{y} }
  \item{qu}{An alternative to \code{th}, specifying the quantile of the \code{y}}
  \item{phi}{ Formula for the log of the scale parameter. Defaults to
      \code{phi = ~ 1} - i.e. no covariates. }
  \item{xi}{ Formula for the shape parameter. Defaults to
      \code{phi = ~ 1} - i.e. no covariates. }
  \item{penalty}{ How to penalize the likelhood. Currently, either
     ``none'', ``gaussian'' or ``lasso'' are the only allowed arguments.
     See the information on the \code{prior} argument, below.}
  \item{prior}{If \code{method = "optimize"}, just an alternative way of
     specifying the pentalty, and only one or neither of \code{penalty}
     and \code{prior} should be given. If \code{method = "simulate"},
     prior must be ``gaussian'' because no other prior distributions
     have been implemented.}
  \item{method}{Should be either ``optimize'' (the default) or ``simulate''.
     The first letter or various abbreviations will do. If ``optimize'' is
     used, the (penalized) likelihood is directly optimized using \code{optim}
     and point estimates (either ML or MAP estimates) are returned with other
     information. If ``simulate'', a Metropolis algorithm is used to simulate
     from the joint posterior distribution of the parameters.}
  \item{start}{Starting values for the parameters, to be passed to \code{optim}.
           If not provided, an exponential distribution is assumed as the starting
           point.}
  \item{priorParameters}{  A list with two components. The first should be 
                a vector of means, the second should be a covariance matrix if the 
                penalty/prior is "gaussian" and a precision matrix is the penalty is "lasso".
                These represent the parameters in the Gaussian prior distribution.
                If not supplied: all default prior means are zero;
                all default prior variances are $10^4$;all covariances
                are zero.  }
  \item{maxit}{ The number of iterations allowed in \code{optim} }
  \item{trace}{ Whether or not to print progress to screen. If \code{method = "optimize"},
       the argument is passed into \code{optim} -- see the help for that
       function. If \code{method = "simulate"}, the argument determines at
       how many steps of the Markov chain the function should tell the user, and
       in this case it defaults to \code{trace = 1000}. }
  \item{iter}{Number of simulations to generate under \code{method = "simulate"}.}
  \item{burn}{ The number of initial steps to be discarded. }
  \item{thin}{ The degree of thinning of the resulting Markov chains. Defaults
               to 1 (no thinning). A value of 0.5 (for example) would result in
               every other step being discarded. }
  \item{jump.const}{ Control parameter for the Metropolis algorithm. }
  \item{verbose}{Whether or not to print progress to screen. Defaults to
                \code{verbose=TRUE}.}
  \item{x, object}{Object of class \code{gpd}, \code{bgpd}, \code{summary.gpd}
                   or \code{summary.bgpd} returned by \code{gpd} or \code{summary.gpd}.}
  \item{digits}{Number of digits for printing.}
  \item{main}{Titles for plots. Should be a vector of length 4.}
  \item{xlab}{In \code{plot.gpd}, x-axis labels for plots. Should be a vector of length 4.}
  \item{nsim}{The number of simulations to use to produce the simulation
              envelope. Defaults to \code{nsim = 1000}}.
  \item{alpha}{A (1 - alpha)\% simulation envelope is produced. Defaults
               to \code{alpha = 0.05}}
  \item{k}{Constant used in calculation of AIC=-2*loglik + k*p, defaults to \code{k=2}.}
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
  \item{dig}{Number of digits to print.}
  \item{...}{Further arguments to be passed to methods.}
}
\details{
  Working with the log of the scale parameter was found to improve the
  stability of computations, and it makes a quadratic penalty more appropriate.
  A quadratic penalty can be thought of as a Gaussian prior distribution, whence
  the terminology of the function.
  
  Some of the internal code is based on the \code{gpd.fit} function
  in the \code{ismev} package and is due to Stuart Coles.

  When a summary or plot is performed, a simulation envelope is
  produced the data, based on quantiles of the fitted
  model. This represents a pointwise (1 - alpha)\%   simulated confidence interval.
  Since the ordered observations will be correlated, if any observation is outside
  the envelope, it is likely that a chain of observations will be outside the
  envelope. Therefore, if the number outside the envelope is a little more than
  alpha\%, that does not immediately imply a serious shortcoming of the fitted model.

  When \code{method = "optimize"}, the \code{plot} function produces diagnostic plots for
  the fitted generalized Pareto model. A PP-plot, QQ-plot,
  histogram with superimposed generalized Pareto density estimate,
  and a return level plot with confidence interval are produced.
  
  The PP-plot and QQ-plot contain simulated pointwise confidence regions.
  The region is a (1 - alpha)\% region based on \code{nsim} simulated
  samples.

  If \code{start} is not provided, the maximum penalized likelhood point
  estimates are computed and used.

  If \code{method = "simulate"}, the simulation is done by a Metropolis
  algorithm.

  When plotting the object, if the chains have converged on the posterior distributions, the trace
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
    If \code{method = "optimize"}, an object of class \code{gpd}.

    \item{convergence}{Output from \code{optim} relating to whether or
                       not the optimizer converged.}
    \item{message}{A message telling the user whether or not convergence
                   was achieved.}
    \item{cov}{The estimated covariance of the parameters in the model.}
    \item{se}{The estimated standard errors of the parameters in the model.}
    \item{threshold}{The threshold of the data above which the GPD model was fit.}
    \item{penalty}{The type of penalty function used, if any.}
    \item{coefficients}{The parameter estimates as computed under maximum
               likelihood or maximum penalized likelihood.}
    \item{rate}{The proportion of observations above the threshold.}
    \item{call}{The call to \code{gpd} that produced the object.}
    \item{y}{The response data.}
    \item{X.phi}{The design matrix for the log of the scale parameter.}
    \item{X.xi}{The design matrix for the shape parameter.}
    \item{loglik}{The value of the optimized log-likelihood.}

    If \code{method = "simulate"}, an object of class \code{bgpd}.
    \item{call}{The call to \code{gpd} that produced the object.}
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

\author{Janet E Heffernan, Harry Southworth}

\examples{
  x <- rnorm(1000)
  mod <- gpd(x, qu = 0.7)
  mod
  par(mfrow=c(2, 2))
  plot(mod)
# Following lines commented out to keep CRAN robots happy
#  mod <- gpd(x, qu=.7, method="sim")
#  mod
#  par(mfrow=c(3, 2))
#  plot(mod)
}
\keyword{ models }

