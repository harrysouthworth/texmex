qgev <- function(p, mu, sigma, xi, lower.tail=TRUE, log.p=FALSE){

    if (log.p == FALSE && (any(p <= 0) || any(p >= 1))) {
        stop("p must lie between 0 and 1 if log.p=FALSE")
    }

    n <- max(length(p), length(sigma), length(xi), length(mu))
    p <- rep(p, length = n)
    sigma <- rep(sigma, length = n)
    xi <- rep(xi, length = n)
    mu <- rep(mu, length = n)

    tol <- 0.0001
    res <- mu - sigma/xi * (1 - (-log(p))^(-xi))
    res[abs(xi) < tol] <- mu - sigma*(log(-log(p)))
    res
}
