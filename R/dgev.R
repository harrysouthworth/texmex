dgev <- function(x, mu, sigma, xi, log.d=FALSE){
    n <- length(x)

    mu <- rep(mu, length=n)
    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)

    tx <- (1 + xi/sigma * (x - mu))

    ld <- -log(sigma) - (1 + 1/xi)*log(tx) - tx^(-1/xi)

    if (!log.d){ exp(ld) }
    else { ld }
}
