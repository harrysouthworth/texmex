coef.bgpd <- function(object, ...){
    apply( object$param, 2, mean )
}
