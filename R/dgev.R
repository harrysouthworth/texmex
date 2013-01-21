dgev <- function(x, sigma, xi, log.d=FALSE){
    n <- length(x)

    mu <- rep(mu, length=n)
    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)

    tx <- (1 + xi/sigma * (x - mu))^(-1/xi)

    if (log.d) { -tx }
    else { exp(-tx) }
}
