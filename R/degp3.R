degp3 <- function(x, kappa=1, sigma, xi, u=0, log.d=FALSE){

  below <- kappa <= 0 | x < u
  x <- pmax((x - u) / sigma, 0)

  xi.ne.0 <- function(x, k, xi, sigma){
      y <- 1 + x * xi
      log(k/sigma) + (k - 1) * log(1 - y^(-1/xi)) - (1/xi + 1) * log(y)
  }
  xi.eq.0 <- function(x, k, xi, sigma){
      log(k/sigma) + (k - 1) * log(1 - exp(-x)) - x
  }

  # Let NaNs be NaNs. Suppress warnings
  res <- suppressWarnings(ifelse(xi == 0, xi.eq.0(x, kappa, xi, sigma), xi.ne.0(x, kappa, xi, sigma)))

  res[below] <- -Inf

  if (log.d) res
  else exp(res)
}
