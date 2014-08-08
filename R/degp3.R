degp3 <- function(x, kappa=1, sigma, xi, u=0, log.d=FALSE){
    res <- log(kappa) + dgpd(x, sigma, xi, u=u, log.d=TRUE) +
        (kappa - 1) * pgpd(x, sigma, xi, u=u, lower.tail=TRUE, log.p=TRUE)

    if (log.d) res
    else exp(res)
}
