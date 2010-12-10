mex <- function(data, which, mth, mqu, dth, dqu,
                penalty="gaussian", maxit=10000,
                trace=0, verbose=FALSE, priorParameters=NULL){

    res1 <- migpd(data=data, mth=mth, mqu=mqu, penalty=penalty,
                  maxit=maxit, trace=trace, verbose=verbose,
                  priorParameters=priorParameters)
    res2 <- mexDependence(x= res1, which=which, dth=dth, dqu=dqu)
    
    res <- list(margins=res1, dependence=res2)
    oldClass(res) <- "mex"
    res
}


print.mex <- function(x, ...){
    cat("Marginal models:\n\n")
    summary(x[[1]])
    cat("Dependence model:\n\n")
    print(x[[2]])
    invisible()
}
show.mex <- print.mex
summary.mex <- function(object, ...){
    print(object, ...)
    invisible(coef(object))
}

plot.mex <- function(x, ...){
    plot(x[[2]])
    invisible()
}

coefficients.mex <- function(object, ...){
    res1 <- coef(object[[1]])
    res2 <- coef(object[[2]])
    list(margins=res1, dependence=res2)
}

coef.mex <- coefficients.mex


plot.mex <- function(x,quantiles=seq(0.1,by=0.2,len=5),col="grey",...){
   if ( class( x ) != "mex" )
		 stop( "you need to use an object with class 'mex'" )

   plot(x$dependence,quantiles=quantiles,col=col,...)
}
 
