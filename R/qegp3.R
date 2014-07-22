qegp3 <- function(p, kappa=1, sigma, xi, u=0, lower.tail=TRUE, log.p=FALSE){
    p <- p^(1/kappa)
    qgpd(p, sigma, xi, u, lower.tail=lower.tail, log.p=log.p)
}
