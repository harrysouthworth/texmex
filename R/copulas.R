edf <- function(x, na.last=NA){
    res <- rank(x) / (length(x) + 1)
    oldClass(res) <- "edf"
    invisible(res)
}

copula <- function(x, na.last=NA){
    stopifnot(is.numeric(x))  
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

test(copula) <- function(){
  fun <- function(d) apply(d,2,function(x)(1:n)[rank(x)])/(1+n)
  n <- 200

  u2 <- cbind(sample(n),sample(n))
  d2 <- fun(u2)

  u3 <- cbind(sample(n),sample(n),sample(n))
  d3 <- fun(u3)
  
  checkEqualsNumeric(d2,copula(u2))
  checkEqualsNumeric(d3,copula(u3))
  checkException(copula(TRUE))
  checkException(copula("text"))
}