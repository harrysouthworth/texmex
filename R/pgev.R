pgev <- function(q, mu, sigma, xi, lower.tail=TRUE, log.p=FALSE){
    tx <- (1 + xi/sigma * (q - mu))^(-1/xi)
    
    tol <- 0.0001
    tx[abs(xi) < tol] <- exp(-(q-mu)/sigma)
    
    exp(-tx)
}
