qegp3 <- function(p, kappa=1, sigma=1, xi, u=0, lower.tail=TRUE){
#    n <- max(length(p), length(k), length(u), length(sigma), length(xi))
#    k <- rep(k, n)
#    p <- rep(p, n)
    
    p <- p^(1/kappa)
    qgpd(p, sigma, xi, u, lower.tail=lower.tail)
}
