#' Information Criteria
#'
#' Compute AIC and (approximate) DIC for \code{evmOpt} objects
#'
#' @aliases AIC.evmSim
#' @param object fit model object
#' @param penalized whether to use the penalized log-likelihood
#' @param nsamp Number of approximate Gaussian sample to use in computing DIC.
#'   Defaults to \code{nsamp=1e3}. Only used when the object has class 'evmOpt'.
#' @param ... other arguments currently ignored
#' @param k numeric, the penalty per parameter to be used; the
#'     default \code{k = 2} is the classical AIC.
#' @details If the object has class 'evmOpt', \code{nsamp} random draws are
#'   made from the Gaussian distribution with mean and covariance inferred from
#'   the model object. The result will be an approximate DIC. Note that AIC should
#'   not be trusted if priors are not flat. For example, if you use a regularizing
#'   prior on xi, say xi ~ N(0, 0.25), AIC can be misleading and DIC should be
#'   preferred. If the object has class 'evmSim', the actual posterior draws are
#'   used in the computation. Also note that sometimes the optimizer returns
#'   an approximatae covariance that is not postive-semidefinite, in which case
#'   the DIC will be reported as NA.
#' @return The AIC and DIC
#' @seealso \code{\link[stats]{AIC}}
#' @importFrom stats AIC logLik
#' @export
AIC.evmOpt <- function(object, penalized=FALSE, nsamp=1e3, ..., k=2){
  ll <- unclass(logLik(object, penalized=penalized))
  aic <- -2*ll + k * attr(ll, 'df')
  attr(aic, "df") <- NULL

  # Get approximate DIC by sampling from Gaussian
  if (all(eigen(object$cov)$values > 0)){ # PSD
    samp <- try(rmvnorm(nsamp, coef(object), object$cov), silent=TRUE) # can throw error if very short tailed so cov mat is singular
    dic <- try(DIC.evm(object, samp), silent=TRUE)

    if (class(dic) == "try-error"){
      dic <- NA
    }
  } else {
    dic <- NA
  }

  c(AIC=aic, DIC=dic)
}

DIC.evm <- function(object, samp){
  ll <- object$family$log.lik(object$data, object$th)

  dev <- -2 * apply(samp, 1, ll)
  dev <- dev[dev < Inf]

  mean(dev) + (mean(dev) + 2 * ll(colMeans(samp)))
}

#' @export
AIC.evmSim <- function(object, ..., k=2){
  samp <- object$param

  res <- DIC.evm(object$map, samp)
  names(res) <- "DIC"
  res
}

#' Log-likelihood for evmOpt objects
#'
#' Return the log-likelihood or penalized log-likelihood for
#' \code{evmOpt} objects.
#'
#' @param object fit model object
#' @param penalized whether to return the penalized log-likelihood
#' @param ... some methods need more arguments
#' @return an object of class \code{logLik}
#' @seealso \code{\link[stats]{logLik}}
#' @importFrom stats logLik
#' @export
logLik.evmOpt <- function(object, penalized=FALSE, ...){
  ll <- if (penalized) {object$ploglik} else {object$loglik}
  attr(ll, 'df') <- length(coef(object))
  class(ll) <- 'logLik'
  ll
}
