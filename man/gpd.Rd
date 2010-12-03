\name{gpd}
\alias{gpd}
\alias{print.gpd}
\alias{summary.gpd}
\alias{plot.gpd}
\alias{coef.gpd}
\alias{AIC.gpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit a generlized Pareto distribution by maximum
    (penalized) likelihood }
\description{
  Fit a generlized Pareto distribution by maximum
    (penalized) likelihood, possibly with explanatory variables.
}
\usage{
gpd(y, data, th, qu, phi = ~1, xi = ~1, penalty = "gaussian",
    start = NULL, priorParameters = NULL, maxit = 10000, trace = 0)
\method{print}{gpd}(x, ...)
\method{summary}{gpd}(object, nsim=1000, alpha=0.05, ...)
\method{coef}{gpd}(object, ...)
\method{plot}{gpd}(x, main=rep(NULL, 4), xlab=rep(NULL, 4), nsim=1000, alpha=0.05, ...)
\method{AIC}{gpd}(object, ..., k=2)
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
     ``none'', ``gaussian'' or ``lasso'' are the only allowed arguments. }
  \item{start}{Starting values for the parameters, to be passed to \code{optim}.
           If not provided, an exponential distribution is assumed as the starting
           point.}
  \item{priorParameters}{  A list with two components. The first should be 
                a vector of means, the second should be a covariance matrix.
                These represent the parameters in the Gaussian prior distribution.
                If not supplied: all default prior means are zero;
                all default prior variances are $10^4$;all covariances
                are zero.  }
  \item{maxit}{ The number of iterations allowed in \code{optim} }
  \item{trace}{ Whether or not to print progress to screen. }
  \item{main}{Titles for plots. Should be a vector of length 4.}
  \itme{xlab}{X-axis labels for plots. Should be a vector of length 4.}
  \item{nsim}{The number of simulations to use to produce the simulation
              envelope. Defaults to \code{nsim = 1000}}.
  \item{alpha}{A (1 - alpha)\% simulation envelope is produced. Defaults
               to \code{alpha = 0.05}}
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

  The \code{plot} function produces diagnostic plots for
  the fitted generalized Pareto model. A PP-plot, QQ-plot,
  histogram with superimposed generalized Pareto density estimate,
  and a return level plot with confidence interval are produced.
  
  The PP-plot and QQ-plot contain simulated pointwise confidence regions.
  The region is a (1 - alpha)\% region based on \code{nsim} simulated
  samples.

}
\value{
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
}

\author{Janet E Heffernan, Harry Southworth}

\examples{
  x <- rnorm(1000)
  mod <- gpd(x, qu = 0.7)
  mod
}
\keyword{ models }

