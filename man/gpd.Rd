\name{gpd}
\alias{gpd}
\alias{gpd.default}
\alias{print.gpd}
\alias{show.gpd}
\alias{summary.gpd}
\alias{plot.gpd}
\alias{coef.gpd}
\alias{AIC.gpd}
\alias{print.bgpd}
\alias{plot.bgpd}
\title{ Generalized Pareto distribution modelling }
\description{Likelihood based modelling and inference for the generalized
Pareto distribution, possibly with explanatory variables.
}
\usage{
gpd(y, data, ...)

\method{gpd}{default}(y, data, th, qu, phi = ~1, xi = ~1, penalty = "gaussian",
    prior = "gaussian", method = "optimize", cov="observed", start = NULL, 
    priorParameters = NULL, maxit = 10000, trace = NULL,
    iter = 40500, burn = 500, thin = 4,
    proposal.dist = c("gaussian", "cauchy"),
    jump.cov, jump.const, verbose = TRUE,...)

\method{print}{gpd}(x, digits=max(3, getOption("digits") - 3), ...)
\method{summary}{gpd}(object, nsim=1000, alpha=0.05, ...)
\method{show}{gpd}(x, digits=max(3, getOption("digits") - 3), ...)

\method{plot}{gpd}(x, main=rep(NULL, 4), xlab=rep(NULL, 4), nsim=1000, alpha=0.05, ...)

\method{AIC}{gpd}(object, ..., k=2)

\method{print}{bgpd}(x, print.seed=FALSE, ...)
\method{plot}{bgpd}(x, which.plots=1:3, density.adjust=2, print.seed=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Either a numeric vector or the name of a variable in \code{data}. }
  \item{data}{ A data frame containing \code{y} and any covariates. }
  \item{th}{ The threshold for \code{y}, exceedances above which will be used to fit the GPD upper tail model. }
  \item{qu}{An alternative to \code{th}, a probability defined such that \code{quantile(y,qu)} equals \code{th}.}
  \item{phi}{ Formula for the log of the scale parameter. Defaults to
      \code{phi = ~ 1} - i.e. no covariates. }
  \item{xi}{ Formula for the shape parameter. Defaults to
      \code{xi = ~ 1} - i.e. no covariates. }
  \item{penalty}{ How to penalize the likelhood. Currently, either
     ``none'', ``gaussian'' or ``lasso'' are the only allowed arguments. 
     If \code{penalty} is "gaussian" or "lasso" then the parameters for 
     the penalization are specified through the \code{prior} argument, see below.}
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
  \item{cov}{How to compute the covariance matrix of the parameters. Defaults to
        \code{cov = "observed"} in which case the observed information matrix
        is used, as given in Appendix A of Davison and Smith. The only other
        option is \code{cov = "numeric"} in which case a numerical approximation
        of the Hessian is used (see the help for \code{optim}). In some cases,
        particularly with small samples, the numerical approximation can be
        quite different from the closed form (\code{cov="observed"}) result, and the value derived from 
        the observed information should be preferred. However, in either case, since the 
        underlying log-likelihood may be far from quadratic for small samples, 
        the resulting estimates of standard errors are liable to approximate 
        poorly the true standard errors. Also see the comments in
        the Details section, below.}
  \item{start}{Starting values for the parameters, to be passed to \code{optim}.
           If not provided, an exponential distribution (shape = 0) is assumed as the starting
           point.}
  \item{priorParameters}{  A list with two components. The first should be 
                a vector of means, the second should be a covariance matrix if the 
                penalty/prior is "gaussian" or "quadratic" and a diagonal precision matrix 
                if the penalty/prior is "lasso", "L1" or "Laplace".
                If \code{method = "simulate"} then these represent the parameters in the Gaussian prior distribution.  If \code{method = "optimize"} then these represent the parameters in the penalty function.
                If not supplied: all default prior means are zero;
                all default prior variances are \eqn{10^4}; all covariances
                are zero.  }
  \item{maxit}{ The number of iterations allowed in \code{optim}. }
  \item{trace}{ Whether or not to print progress to screen. If \code{method = "optimize"},
       the argument is passed into \code{optim} -- see the help for that
       function. If \code{method = "simulate"}, the argument determines at
       how many steps of the Markov chain the function should tell the user, and
       in this case it defaults to \code{trace = 1000}. }
  \item{iter}{Number of simulations to generate under \code{method = "simulate"}.
              Defaults to 40000.}
  \item{burn}{ The number of initial steps to be discarded. }
  \item{thin}{ The degree of thinning of the resulting Markov chains. Defaults
               to 4 (one in every 4 steps is retained).}
  \item{proposal.dist}{The proposal distribution to use, either
    multivariate gaussian or a multivariate Cauchy.}
  \item{jump.cov}{Covariance matrix for proposal distribution of Metropolis algorithm.  This is scaled by \code{jump.const}.}
  \item{jump.const}{ Control parameter for the Metropolis algorithm. }
  \item{verbose}{Whether or not to print progress to screen. Defaults to
                \code{verbose=TRUE}.}
 \item{x, object}{Object of class \code{gpd}, \code{bgpd}, \code{summary.gpd}
                   or \code{summary.bgpd} returned by \code{gpd} or \code{summary.gpd}.}
  \item{digits}{Number of digits for printing.}
  \item{main}{In \code{plot} method for class \code{gpd}, titles for diagnostic plots. Should be a vector of length 4, with values corresponding to the character strings to appear on the titles of the pp- qq- return level and density estimate plots respectively.}
  \item{xlab}{As for \code{main} but labels for x-axes rather than titles.}
  \item{nsim}{In \code{plot} and \code{summary} methods for class \code{gpd}. The number of replicates to be simulated to produce the simulated tolerance intervals. Defaults to \code{nsim = 1000}}.
  \item{alpha}{In \code{plot} and \code{summary} methods for class \code{gpd}. A (1 - alpha)\% simulation envelope is produced. Defaults to \code{alpha = 0.05}}
  \item{k}{Constant used in calculation of AIC=-2*loglik + k*p, defaults to \code{k=2}.}
  \item{print.seed}{Whether or not to print the seed used in the simulations,
                    or to annotate the plots with it. Defaults to
                    \code{print.seed=FALSE}.}
  \item{which.plots}{In \code{plot} method for class \code{bgpd}. Which plots to produce. Option 1 gives kernel density estimates,
                     2 gives traces of the Markov chains with superimposed
                     cumulative means, 3 gives autocorrelation functions.
                     Defaults to \code{which.plots=1:3}.}
  \item{density.adjust}{In \code{plot} method for class \code{bgpd}. Passed into \code{density}. Controls the amount of
                     smoothing of the kernel density estimate. Defaults to
                     \code{density.adjust=2}.}
  \item{...}{Further arguments to be passed to methods.}
}
\details{
  We use the following parameterisation of the GPD:
  
  \deqn{P(X \le x) = 1 - \left(1 + \frac{\xi x}{\sigma}\right)^{-1/\xi}}{%
        P(X \le x) = 1 - (1 + \xi x / \sigma)^(-1/\xi)}
  for \eqn{x \ge 0} and \eqn{1 + \xi x / \sigma \ge 0.}  The scale parameter 
  is sigma (\eqn{\sigma}) and the shape parameter is xi (\eqn{\xi}).
  
  Working with the log of the scale parameter improves the
  stability of computations, making  a quadratic penalty more 
  appropriate and enabling the inclusion of covariates in the model for the scale parameter, which must remain positive.  We therefore work with \eqn{\phi}=log(\eqn{\sigma}).  All specification of priors or penalty functions refer to \eqn{\phi} rather than \eqn{\sigma}.  A quadratic penalty can be 
  thought of as a Gaussian prior distribution, whence the terminology of the function.
  
  Parameters of the GPD fitted to excesses above threshold \code{th} are estimated by using penalized maximum likelihood (\code{method="optimize"}), or by simulating from the posterior distribution of the model parameters using a Metropolis algorithm (\code{method="simulate"}).  In the latter case, \code{start} is used as a starting value for the Metropolis algorithm; in its absence, the maximum penalized likelhood point estimates are computed and used.

  When a summary or plot is performed, a pointwise (1 - alpha)\% tolerance envelope is simulated, based on quantiles of the fitted
  model.  Since the ordered observations will be correlated, if any observation is outside
  the envelope, it is likely that a chain of observations will be outside the
  envelope. Therefore, if the number outside the envelope is a little more than
  alpha\%, that does not immediately imply a serious shortcoming of the fitted model.

  When \code{method = "optimize"}, the \code{plot} function produces 
  diagnostic plots for the fitted generalized Pareto model. These differ 
  depending on whether or not there are covariates in the model.  If there 
  are no covariates then the diagnostic plots are PP- and QQ-plots, a 
  return level plot (produced by \code{plotrl.gpd}) and a histogram of the data with superimposed generalized 
  Pareto density estimate.  These are all calculated using the data on 
  the original scale. If there are covariates in the model then the 
  diagnostics consist of PP- and QQ- plots calculated by using the model 
  residuals (which will be standard exponential devaiates under the GPD model), 
  and plots of residuals versus fitted model parameters.
  
  The PP- and QQ-plots show simulated pointwise tolerance intervals.
  The region is a (1 - alpha)\% region based on \code{nsim} simulated
  samples.

  When \code{method = "simulate"} the \code{plot} function produces 
  diagnostic plots for the Markov chains used to simulate from the posterior 
  distributions for the model parameters.  If the chains have converged on 
  the posterior distributions, the trace plots should look like "fat hairy 
  caterpillars" and their cumulative means should converge rapidly. Moreover, 
  the autocorrelation functions should converge quickly to zero.

  When \code{method = "simulate"} the \code{print} and \code{summary} functions 
  give posterior means and standard deviations. Posterior means are also 
  returned by the \code{coef} method. Depending on what you want to do 
  and what the posterior distributions look like (use \code{plot} method)
  you might want to work with quantiles of the posterior distributions 
  instead of relying on standard errors.
  
 
}
\value{
    If \code{method = "optimize"}, an object of class \code{gpd}:

    \item{convergence}{Output from \code{optim} relating to whether or
                       not the optimizer converged.}
    \item{message}{A message telling the user whether or not convergence
                   was achieved.}
    \item{threshold}{The threshold of the data above which the GPD model was fit.}
    \item{penalty}{The type of penalty function used, if any.}
    \item{coefficients}{The parameter estimates as computed under maximum
                        likelihood or maximum penalized likelihood.}
    \item{rate}{The proportion of observations above the threshold.}
    \item{call}{The call to \code{gpd} that produced the object.}
    \item{y}{The response data above the threshold for fitting.}
    \item{X.phi}{The design matrix for the log of the scale parameter.}
    \item{X.xi}{The design matrix for the shape parameter.}
    \item{priorParameters}{See above.}
    \item{data}{The original data (above and below the threshold for fitting).}
    \item{residuals}{Data above the threshold for fitting after transformation 
                     to standard exponential scale by using the fitted GPD.}
    \item{loglik}{The value of the optimized log-likelihood.}
    \item{cov}{The estimated covariance of the parameters in the model.}
    \item{se}{The estimated standard errors of the parameters in the model.}
    
    

    If \code{method = "simulate"}, an object of class \code{bgpd}:

    \item{call}{The call to \code{gpd} that produced the object.}
    \item{threshold}{The threshold above which the model was fit.}
    \item{map}{The point estimates found by maximum penalized likelihood
               and which were used as the starting point for the Markov
               chain.  This is of class \code{gpd} and methods for this 
               class (such as resid and plot) may be useful.}
    \item{burn}{The number of steps of the Markov chain that are to be 
                treated as the burn-in and not used in inferences.}
    \item{thin}{The degree of thinning used.}
    \item{chains}{The entire Markov chain generated by the Metropolis
                  algorithm.}
    \item{y}{The response data above the threshold for fitting.}  
    \item{data}{The original data (above and below the threshold for fitting).}
    \item{X.phi}{The design matrix for the log of the scale parameter.}
    \item{X.xi}{The design matrix for the log of the shape parameter.}
    \item{acceptance}{The proportion of proposals that were accepted by
                      the Metropolis algorithm.}
    \item{seed}{The seed used by the random number generator.}
    \item{param}{The remainder of the chain after deleting the burn-in
                 and applying any thinning.}
                 
    There are summary, plot, print and coefficients methods available for these classes.

}
\note{   When there are estimated values of \code{xi <= -0.5}, the regularity
  conditions of the likelihood break down and inference based on approximate
  standard errors cannot be performed. In this case, the most fruitful
  approach to inference appears to be by the bootstrap. It might be possible
  to simulate from the posterior, but finding a good proposal distribution
  might be difficult and you should take care to get an acceptance rate
  that is reasonably high (around 40\% when there are no covariates, lower
  otherwise).
  }

\seealso{ \code{\link{rl.gpd}}, \code{\link{predict.gpd}}, \code{\link{gpd.declustered}}.}

\author{Janet E. Heffernan, Harry Southworth.   Some of the internal code is based on the \code{gpd.fit} function
  in the \code{ismev} package and is due to Stuart Coles.}

\references{A. C. Davison and R. L. Smith, Models for exceedances over
            high thresholds, Journal of the Royal Statistical Society B,
            53, 393 -- 442, 1990}

\examples{
  x <- rnorm(1000)
  mod <- gpd(x, qu = 0.7)
  mod
  par(mfrow=c(2, 2))
  plot(mod)

  x <- runif(100,-0.2,0.2)
  data <- data.frame(x=x,y=rgpd(100,sigma=exp(3 + 2*x),xi=x))
  mod <- gpd(y, data, phi = ~x, xi = ~x, th = 0)
  plot(mod)

# Following lines commented out to keep CRAN happy
#  mod <- gpd(x, qu=.7, method="sim")
#  mod
#  par(mfrow=c(3, 2))
#  plot(mod)
}
\keyword{ models }

