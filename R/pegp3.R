pegp3 <- function(q, kappa=1, sigma, xi, u=0, lower.tail=TRUE, log.p=FALSE){
    res <- pgpd(q, sigma=sigma, xi=xi, u=u, lower.tail=lower.tail, log.p=log.p)
    res^kappa
}
