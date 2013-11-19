\name{evm}
\alias{evm}
\alias{evm.default}
\alias{evmSim}
\alias{texmexMetropolis}
\alias{texmexCheckMap}
\alias{texmexJumpConst}
\alias{initRNG}
\alias{evmBoot}
\alias{print.evmOpt}
\alias{summary.evmOpt}
\alias{plot.evmOpt}
\alias{coef.evmOpt}
\alias{AIC.evmOpt}
\alias{print.evmSim}
\alias{summary.evmSim}
\alias{plot.evmSim}
\alias{print.evmBoot}
\alias{summary.evmBoot}
\alias{plot.evmBoot}
\title{Extreme value modelling}
\description{Likelihood based modelling and inference for extreme value models,
possibly with explanatory variables.}
\usage{
evm(y, data, family, ...)

\method{evm}{default}(y, data, family=gpd, th=-Inf, qu,
    ..., penalty = NULL,
    prior = "gaussian", method = "optimize", cov="observed", start = NULL,
    priorParameters = NULL, maxit = 10000, trace = NULL,
    iter = 40500, burn = 500, thin = 4,
    proposal.dist = c("gaussian", "cauchy"),
    jump.cov, jump.const=NULL, R=100, verbose = TRUE)

\method{print}{evmOpt}(x, digits=max(3, getOption("digits") - 3), ...)
\method{summary}{evmOpt}(object, nsim=1000, alpha=0.05, ...)
\method{plot}{evmOpt}(x, main=rep(NULL, 4), xlab=rep(NULL, 4),
           nsim=1000, alpha=0.05, ...)
\method{AIC}{evmOpt}(object, penalized=FALSE, ..., k=2)

\method{print}{evmSim}(x, print.seed=FALSE, ...)
\method{plot}{evmSim}(x, which.plots=1:3, density.adjust=2, print.seed=FALSE, ...)
texmexMetropolis(x, log.lik, proposals, verbose, trace)
texmexJumpConst(j, map)
}

\arguments{
  \item{y}{ Either a numeric vector or the name of a variable in \code{data}. }
  \item{data}{ A data frame containing \code{y} and any covariates. }
  \item{family}{An object of class 'texmexFamily'. Defaults to
           \code{family=gpd} and a generalized Pareto distribution is
           fit to the data. Alternatively the family could be \code{gev},
           resulting in a generalized extreme value distribution being
           fit. No other families are currently available in texmex,
           but users may write their own.}
  \item{th}{ For threshold excess models (such as when
           \code{family=gpd}), the threshold for \code{y},
           exceedances above which will be used to fit the upper tail
           model. Note that if you have already thresholded your data
           and want to model all of \code{y}, you still need to specify
           \code{th}.}
  \item{qu}{An alternative to \code{th}, a probability defined
         such that \code{quantile(y,qu)} equals \code{th}.}
  \item{...}{ In \code{evm}, formulae for the parameters in the
           family, e.g. \code{phi = ~ x}. If none are specified,
           they all default to \code{~1}.}
  \item{penalty}{ How to penalize the likelhood. Currently, either
           "none"", "gaussian"" or "lasso"" are the only allowed
           values. If \code{penalty} is "gaussian" or "lasso" then
           the parameters for the penalization are specified through
           the \code{priorParameters} argument. See below. Defaults to
           \code{penalty=NULL} and applies maximum likelihood estimation.}
  \item{prior}{If \code{method = "optimize"}, just an alternative way of
           specifying the pentalty, and only one or neither of \code{penalty}
           and \code{prior} should be given. If \code{method = "simulate"},
           prior must be "gaussian" because no other prior distributions
           have been implemented.}
  \item{method}{Should be either "optimize" (the default), "simulate"
           or "bootstrap".
           The first letter or various abbreviations will do. If 'optimize' is
           used, the (penalized) likelihood is directly optimized using \code{optim}
           and point estimates (either ML or MAP estimates) are returned with other
           information. If "simulate", a Metropolis algorithm is used to simulate
           from the joint posterior distribution of the parameters. If "bootstrap",
           a parametric boostrap is performed.}
  \item{cov}{How to compute the covariance matrix of the parameters. Defaults to
           \code{cov = "observed"} in which case the observed information matrix
           is used, if the \code{info} element of the \code{texmexFamily} object
           is present. Note that currently, this is not implemented for \code{gev}. The only other
           option is \code{cov = "numeric"} in which case a numerical approximation
           of the Hessian is used (see the help for \code{optim}). In some cases,
           particularly with small samples, the numerical approximation can be
           quite different from the closed form (\code{cov="observed"})
           result, and the value derived from the observed information
           should be preferred. However, in either case, since the
           underlying log-likelihood may be far from quadratic for small samples,
           the resulting estimates of standard errors are liable to approximate
           poorly the true standard errors. Also see the comments in
           the Details section, below.}
  \item{start}{Starting values for the parameters, to be passed to \code{optim}.
           If not provided, the function will use the \code{start} element
           of the \code{texmexFamily} object if it exists.}
  \item{priorParameters}{  A list with two components. The first should be
           a vector of means, the second should be a covariance matrix if the
           penalty/prior is "gaussian" or "quadratic" and a diagonal
           precision matrix
           if the penalty/prior is "lasso", "L1" or "Laplace".
           If \code{method = "simulate"} then these represent the
           parameters in the Gaussian prior distribution.  If
           \code{method = 'optimize'} then these represent the
           parameters in the penalty function.
           If not supplied: all default prior means are zero;
           all default prior variances are \eqn{10^4}; all covariances
           are zero.  }
  \item{maxit}{ The number of iterations allowed in \code{optim}. }
  \item{trace}{ Whether or not to print progress to screen. If
           \code{method = "optimize"},
           the argument is passed into \code{optim} -- see the help for that
           function. If \code{method = "simulate"}, the argument determines at
           how many steps of the Markov chain the function should tell the user, and
           in this case it defaults to \code{trace = 10000}. }
  \item{iter}{Number of simulations to generate under \code{method = "simulate"}.
           Defaults to 40500.}
  \item{burn}{ The number of initial steps to be discarded. Defaults to 500.}
  \item{thin}{ The degree of thinning of the resulting Markov chains. Defaults
           to 4 (one in every 4 steps is retained).}
  \item{proposal.dist}{The proposal distribution to use, either
           multivariate gaussian or a multivariate Cauchy.}
  \item{jump.cov}{Covariance matrix for proposal distribution of
           Metropolis algorithm.  This is scaled by \code{jump.const}.}
  \item{jump.const, j}{ Control parameter for the Metropolis algorithm. }
  \item{verbose}{Whether or not to print progress to screen. Defaults to
           \code{verbose=TRUE}.}
 \item{x, object}{Object of class \code{evmOpt}, \code{evmSim}, \code{summary.evmOpt}
           or \code{summary.evmSim} returned by \code{evmSim} or
           \code{summary.evmSim}. In \code{texmexMetropolis}, \code{x}
           is a matrix, the first row of which is used as the starting
           point of the simulation.}
  \item{digits}{Number of digits for printing.}
  \item{main}{In \code{plot} method for class \code{evmSim}, titles for
           diagnostic plots. Should be a vector of length 4, with values
           corresponding to the character strings to appear on the titles of
           the pp, qq, return level, and density estimate plots respectively.}
  \item{xlab}{As for \code{main} but labels for x-axes rather than titles.}
  \item{nsim}{In \code{plot} and \code{summary} methods for class
           \code{evmSim}. The number of replicates to be simulated to produce
           the simulated tolerance intervals. Defaults to \code{nsim = 1000}.}
  \item{alpha}{In \code{plot} and \code{summary} methods for class \code{evmSim}.
           A 100(1 - alpha)\% simulation envelope is produced. Defaults to
           \code{alpha = 0.05}}
  \item{penalized}{Whether to compute AIC using the penalized log-likelihood or
        the true log-liklihood. Defaults to \code{penalized=FALSE} and uses the
        true log-likelihood.}
  \item{k}{Constant used in calculation of AIC=-2*loglik + k*p,
           defaults to \code{k=2}.}
  \item{print.seed}{Whether or not to print the seed used in the simulations,
           or to annotate the plots with it. Defaults to
           \code{print.seed=FALSE}.}
  \item{which.plots}{In \code{plot} method for class \code{evmSim}. Which plots
           to produce. Option 1 gives kernel density estimates,
           2 gives traces of the Markov chains with superimposed
           cumulative means, 3 gives autocorrelation functions.
           Defaults to \code{which.plots=1:3}.}
  \item{density.adjust}{In \code{plot} method for class \code{evmSim}.
           Passed into \code{density}. Controls the amount of
           smoothing of the kernel density estimate. Defaults to
           \code{density.adjust=2}.}
  \item{log.lik, proposals}{The log-likelihood for the model, and a matrix of random
           numbers drawn from the proposal distribution.}
  \item{map}{In \code{texmexJumpConst}, a model fit by maximum
           (penalized) likelihood. The function is not intended to be called by
           the end user.}
  \item{R}{The number of parametric bootstrap samples to run when
           \code{method = "bootstrap"} is requested. Defaults to 100.}
}
\details{
  In previous versions of \code{texmex}, the main modelling
  function was \code{gpd} and models using the generalized
  Pareto distribution (GPD) were used. In this version, the main
  modelling function is now \code{evm} (extreme value model)
  and the distribution to be used is specified through the
  \code{family} argument.

  Object \code{gpd} is now an object of class \code{texmexFamily}.
  Currently, the only other \code{texmexFamily} object available
  is \code{gev} which results in fitting a generalized extreme
  value (GEV) distribution to the data.

  See Coles (2001) for an introduction to extreme value modelling and
  the GPD and GEV models.

  For the GPD model, we use the following parameterisation of evm:

  \deqn{P(Y \le y) = 1 - \left(1 + \frac{\xi y}{\sigma}\right)^{-1/\xi}}{%
        P(Y \le y) = 1 - (1 + \xi y / \sigma)^(-1/\xi)}
  for \eqn{y \ge 0} and \eqn{1 + \xi y / \sigma \ge 0.}

  For the GEV model, we use:
  
  \deqn{P(Y \le y) = \exp \left\{ -\left[ 1 + \xi \left( \frac{y - \mu}{\sigma}\right) \right] ^{-1/\xi} \right\} }{%
        P(Y \le y) = exp{-[1 + \xi(y - \mu)/\sigma]^{1/\xi}}}

  In each case, the scale parameter is sigma (\eqn{\sigma}) and the
  shape parameter is xi (\eqn{\xi}). The GEV distribution also has
  location parameter mu (\eqn{\mu}).

  Working with the log of the scale parameter improves the
  stability of computations, makes  a quadratic penalty more
  appropriate and enables the inclusion of covariates in the model
  for the scale parameter, which must remain positive.  We therefore
  work with \eqn{\phi}=log(\eqn{\sigma}).  All specification of
  priors or penalty functions refer to \eqn{\phi} rather
  than \eqn{\sigma}.  A quadratic penalty can be thought of as a
  Gaussian prior distribution, whence the terminology of the function.

  Parameters of the evm are estimated by using maximum (penalized)
  likelihood (\code{method = "optimize"}), or by simulating from the
  posterior distribution of the model parameters using a
  Metropolis algorithm (\code{method = "simulate"}).  In the
  latter case, \code{start} is used as a starting value for the
  Metropolis algorithm; in its absence, the maximum penalized likelhood
  point estimates are computed and used.

  A boostrap approach is also available (\code{method = "bootstrap"}).
  This runs a parametric bootstrap, simulating from the model fit
  by optimization.

  When \code{method = "optimize"}, the \code{plot} function produces
  diagnostic plots for the fitted model. These differ
  depending on whether or not there are covariates in the model.  If there
  are no covariates then the diagnostic plots are PP- and QQ-plots, a
  return level plot (produced by \code{plotrl.evmSim}) and a histogram
  of the data with superimposed density estimate.  These are all calculated
  using the data on the original scale. If there are covariates in the
  model then the diagnostics consist of PP- and QQ- plots calculated
  by using the model residuals (which will be standard exponential
  devaiates under the GPD model and standard Gumbel deviates under
  the GEV model), and plots of residuals versus fitted model parameters.

  The PP- and QQ-plots show simulated pointwise tolerance intervals.
  The region is a 100(1 - alpha)\% region based on \code{nsim} simulated
  samples.

  When \code{method = "simulate"} the \code{plot} function produces
  diagnostic plots for the Markov chains used to simulate from the posterior
  distributions for the model parameters. If the chains have converged on
  the posterior distributions, the trace plots should look like "fat hairy
  caterpillars" and their cumulative means should converge rapidly. Moreover,
  the autocorrelation functions should converge quickly to zero.

  When \code{method = "bootstrap"}, the \code{plot} function plots
  the bootstrap distributions of the parameters.

  When \code{method = "simulate"} the \code{print} and \code{summary} functions
  give posterior means and standard deviations. Posterior means are also
  returned by the \code{coef} method. Depending on what you want to do
  and what the posterior distributions look like (use \code{plot} method)
  you might want to work with quantiles of the posterior distributions
  instead of relying on standard errors.

  When \code{method = "bootstrap"}, summaries of the bootstrap distribution
  and the bootstrap estimate of bias are displayed.
}
\value{
    If \code{method = "optimize"}, an object of class \code{evmOpt}:

    \item{call}{The call to \code{evmSim} that produced the object.}
    \item{data}{The original data (above and below the threshold for fitting if
          a distribution for threshold excesses has been used). In detail,
          \code{data} is a list with elements \code{y} and \code{D}. \code{y}
          is the response variable and \code{D} is a list containing the design
          matrices implied by any formlae used in the call to \code{evm}.}
    \item{convergence}{Output from \code{optim} relating to whether or
                       not the optimizer converged.}
    \item{message}{A message telling the user whether or not convergence
                   was achieved.}
    \item{threshold}{The threshold of the data above which the evmSim model was fit.}
    \item{penalty}{The type of penalty function used, if any.}
    \item{coefficients}{The parameter estimates as computed under maximum
                        likelihood or maximum penalized likelihood.}
    \item{rate}{The proportion of observations above the threshold. If the
            model is not a threshold exceedance model (e.g. the GEV
            model), the rate will be 1.}

    \item{priorParameters}{See above.}

    \item{residuals}{Residuals computed using the residual function in the
            \code{texmexFamily} object, if any.}
    \item{ploglik}{The value of the optimized penalized log-likelihood.}
    \item{loglik}{The value of the optimized (unpenalized) log-likelihood. If
          \code{penalty='none'} is used, this will be identical to \code{ploglik}, above.}
    \item{cov}{The estimated covariance of the parameters in the model.}
    \item{se}{The estimated standard errors of the parameters in the model.}
    \item{xlevels}{A named list containing a named list for each design matrix (main parameter)
          in the model. Each list contians an element named after each factor in the
          linear predictor for the respective design matrix. These are used by the \code{predict}
          method to ensure all factor levels are known, even if they don't appear in
          \code{newdata}.}

    If \code{method = "simulate"}, an object of class \code{evmSim}:

    \item{call}{The call to \code{evmSim} that produced the object.}
    \item{threshold}{The threshold above which the model was fit.}
    \item{map}{The point estimates found by maximum penalized likelihood
               and which were used as the starting point for the Markov
               chain.  This is of class \code{evmOpt} and methods for this
               class (such as resid and plot) may be useful.}
    \item{burn}{The number of steps of the Markov chain that are to be
                treated as the burn-in and not used in inferences.}
    \item{thin}{The degree of thinning used.}
    \item{chains}{The entire Markov chain generated by the Metropolis
                  algorithm.}
    \item{y}{The response data above the threshold for fitting.}
    \item{seed}{The seed used by the random number generator.}
    \item{param}{The remainder of the chain after deleting the burn-in
                 and applying any thinning.}

    If \code{method = "bootstrap"}, an object of class \code{evmBoot}:
    
    \item{call}{The call to \code{evmBoot} that produced the object.}
    \item{replicates}{The parameter estimates from the bootstrap fits.}
    \item{map}{The fit by by maximum penalized likelihood to the orginal data.
               This is of class \code{evmOpt} and methods for this
               class (such as resid and plot) may be useful.}

    There are summary, plot, print and coefficients methods available for these classes.

}
\note{   For both GPD and GEV models, when there are estimated values
  of \eqn{\xi \le -0.5}, the regularity
  conditions of the likelihood break down and inference based on approximate
  standard errors cannot be performed. In this case, the most fruitful
  approach to inference appears to be by the bootstrap. It might be possible
  to simulate from the posterior, but finding a good proposal distribution
  might be difficult and you should take care to get an acceptance rate
  that is reasonably high (around 40\% when there are no covariates, lower
  otherwise).
}

  \seealso{ \code{\link{rl.evmOpt}}, \code{\link{predict.evmOpt}},
    \code{\link{evm.declustered}}.}

\author{Janet E. Heffernan, Harry Southworth. Some of the internal code is based on
  the \code{gpd.fit} function in the \code{ismev} package and is due to Stuart Coles.}

\references{S. Coles. An Introduction to Statistical Modelling of Extreme
            Values. Springer, 2001.}

\examples{
  mod <- evm(rain, th=30)
  mod
  par(mfrow=c(2, 2))
  plot(mod)

#  mod <- evm(rain, th=30, method="sim")
#  par(mfrow=c(3, 2))
#  plot(mod)

  mod <- evm(SeaLevel, data=portpirie, family=gev)
  mod
  plot(mod)

#  mod <- evm(SeaLevel, data=portpirie, family=gev, method="sim")
#  par(mfrow=c(3, 3))
#  plot(mod)
}
\keyword{ models }

