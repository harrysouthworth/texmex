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
    print(summary(x[[1]]))
    print(x[[2]])
    invisible()
}
show.mex <- print.mex
summary.mex <- function(object, ...){
    print(summary(object[[1]]))
    print(object[[2]])
    invisible(coef(object))
}

plot.mex <- function(x, ...){
    plot(x[[2]])
    invisible()
}

coefficients.mex <- function(x, ...){
    res1 <- coef(x[[1]])
    res2 <- coef(x[[2]])
    list(margins=res1, dependence=res2)
}
coef.mex <- coefficients.mex
