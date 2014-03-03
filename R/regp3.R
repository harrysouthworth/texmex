regp3 <- function(n, kappa=1, sigma, xi, u=0){
  if (kappa != 1){
    x <- 1 - runif(n)^(1/kappa)
    (x^(-xi) - 1) * sigma / xi + u
  }
  else {
    rgpd(n, sigma, xi, u)
  }
}


