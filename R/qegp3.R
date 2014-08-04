qegp3 <- function(p, kappa=1, sigma, xi, u=0, lower.tail=TRUE, log.p=FALSE){
  if (lower.tail){
    p <- p^(1/kappa)
    qgpd(p, sigma, xi, u, lower.tail=lower.tail, log.p=log.p)
  } else {
    # If we want the upper tail, find the lower tail probability, transform,
    # find the GPD lower tail, then reconvert to upper tail
    p <- (1- p)^(1/kappa)
    res <- qgpd(p, sigma, xi, u, lower.tail=TRUE, log.p=FALSE)
    if (log.p) log(1 - exp(res))
    else 1 - res
  }
}
