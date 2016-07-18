#' Akaike's Information Criterion
#'
#' Compute AIC for \code{evmOpt} objects
#'
#' @param object fit model object
#' @param penalized whether to use the penalized log-likelihood
#' @param ... other arguments currently ignored
#' @param k numeric, the penalty per parameter to be used; the
#'     default \code{k = 2} is the classical AIC.
#' @return the AIC
#' @seealso \code{\link[stats]{AIC}}
#' @importFrom stats AIC logLik
#' @export
AIC.evmOpt <- function(object, penalized=FALSE, ..., k=2){
  ll <- unclass(logLik(object, penalized=penalized))
  aic <- -2*ll + k * attr(ll, 'df')
  attr(aic, "df") <- NULL
  aic
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
