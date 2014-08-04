regp3 <- function(n, kappa=1, sigma, xi, u=0){

  kappa <- rep(kappa, length.out=n)
  sigma <- rep(sigma, length.out=n)
  xi <- rep(xi, length.out=n)
  u <- rep(u, length.out=n)

  nk1 <- sum(kappa != 1)

  x <- 1 - runif(nk1)^(1/kappa)
  
  xi.neq.0 <- function(x, xi, sigma, u) (x^(-xi) - 1) * sigma / xi + u
  xi.eq.0 <- function(x, xi, sigma, u) log(x) * sigma / xi + u
  
  res <- ifelse(xi == 0, xi.eq.0(x, xi, sigma, u), xi.neq.0(x, xi, sigma, u))
  
  if (nk1 < n){
    res2 <- rgpd(n - nk1, sigma, xi, u)
    r <- rep(0, n)
    r[kappa != 1] <- res
    r[kappa == 1] <- res2
    r
  } else res
}


