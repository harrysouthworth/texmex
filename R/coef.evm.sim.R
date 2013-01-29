coef.evm.sim <- function(object, ...){
    res <- apply( object$param, 2, mean )
    names(res) <- c(paste("phi:", colnames(object$X.phi)), 
        paste("xi:", colnames(object$X.xi)))
    res
}
coefficients.evm.sim <- coef.evm.sim

