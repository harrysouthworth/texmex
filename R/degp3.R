degp3 <- function(x, kappa=1, sigma=1, xi, u=0, log.d=FALSE){
  
    below <- kappa <= 0 | x < u
    
    x <- pmax((x - u) / sigma, 0)
    
    n <- max(length(x), length(kappa), length(u), length(sigma), length(xi))
#    sigma <- rep(sigma, n)
    xi <- rep(xi, n)
    kappa <- rep(kappa, n)
#    u <- rep(u, n)

    xi.ne.0 <- function(x, k, xi, sigma){
        y <- 1 + x * xi
        log(k/sigma) + (k - 1) * log(1 - y^(-1/xi)) - (1/xi + 1) * log(y)
    }
    xi.eq.0 <- function(x, k, xi, sigma){
        log(k/sigma) + (k - 1) * log(1 - exp(-x)) - x
    }
#browser()
    # Let NaNs be NaNs. Suppress warnings
    res <- suppressWarnings(ifelse(xi == 0, xi.eq.0(x, kappa, xi, sigma), xi.ne.0(x, kappa, xi, sigma)))

    res[below] <- -Inf
    
    if (log.d){ res }
    else { exp(res) }
}
