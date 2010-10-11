edf <- function(x, na.last=NA){
    res <- rank(x) / (length(x) + 1)
    oldClass(res) <- "edf"
    invisible(res)
}

copula <- function(x, na.last=NA){
    res <- apply(x, 2, edf)
    oldClass(res) <- "copula"
    res
}

plot.copula <- function(x, jitter.=FALSE, ...){
    if (jitter.){
        x <- apply(x, 2, jitter)
    }
    pairs(x, ...)
    invisible()
}

