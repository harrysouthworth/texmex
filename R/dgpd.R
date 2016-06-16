#' Density, cumulative density, quantiles and random number generation for the
#' generalized Pareto distribution
#' 
#' Density, cumulative density, quantiles and random number generation for the
#' generalized Pareto distribution
#' 
#' Random number generation is done by transformation of a standard
#' exponential.
#' 
#' @param x,q,p Value, quantile or probability respectively.
#' @param n Number of random numbers to simulate.
#' @param sigma Scale parameter.
#' @param xi Shape parameter.
#' @param u Threshold
#' @param log.d,log.p Whether or not to work on the log scale.
#' @param lower.tail Whether to return the lower tail.
#' @author Janet E Heffernan, Paul Metcalfe, Harry Southworth
#' @keywords models
#' @examples
#' 
#'   x <- rgpd(1000, sigma=1, xi=.5)
#'   hist(x)
#'   x <- rgpd(1000, sigma=exp(rnorm(1000, 1, .25)), xi=rnorm(1000, .5, .2))
#'   hist(x)
#'   plot(pgpd(x, sigma=1, xi=.5))
#' 
#' @export dgpd
dgpd <- function(x, sigma, xi, u = 0, log.d=FALSE ) {
  ## record the things below bounds for later
  below.bounds <- x < u
  ## shift and scale
  x <- pmax((x - u) / sigma, 0)

  xix <- .specfun.safe.product(xi, x)
  logrel <- .log1prel(xix) * x

  log.density <- -log(sigma) - log1p(xix) - logrel
  log.density[below.bounds] <- -Inf
  log.density[xix==-1] <- -Inf

  if (!log.d) {
    exp(log.density)
  } else {
    log.density
  }
}


