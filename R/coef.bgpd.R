coef.bgpd <- function(object, ...){
    res <- apply( object$param, 2, mean )
    names(res) <- c(paste("phi:", colnames(object$X.phi)), 
        paste("xi:", colnames(object$X.xi)))
    res
}
coefficients.bgpd <- coef.bgpd

