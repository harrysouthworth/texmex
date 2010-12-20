edf <- function(x, na.last=NA){
    res <- rank(x) / (length(x) + 1)
    oldClass(res) <- "edf"
    invisible(res)
}

copula <- 
function (x, na.last = NA) {
    theCall <- match.call()
    
    if (is.data.frame(x)){
        really.numeric <- function(x){
            if (! class(x) %in% c("integer", "numeric")){ FALSE }
            else { TRUE }
        }

        wh <- sapply(x, really.numeric)
    
        if (sum(wh) == 0){
            stop("x contains no numeric columns")
        }
    
        if (sum(wh) < length(wh)){
            warning(paste("Some variables have been dropped:", paste(colnames(x)[!wh], collapse=", ")))
        }

        x <- as.matrix(x[, wh])
    } # Close if
    
    else if (!is.matrix(x)){
        stop("x should be a matrix or a data.frame with some numeric columns")
    }
    
    res <- apply(x, 2, edf)

    res <- list(call=theCall, copula=res)
    oldClass(res) <- "copula"
    res
}


print.copula <- function(x, ...){
    print(x$call)
    cat("A copula of", ncol(x$copula), "variables.\n")
    invisible()
}

show.copula <- print.copula

summary.copula <- function(object, ...){
    print(object$call)
    cat("A copula of", ncol(object$copula), "variables.\n")
    invisible()
}

plot.copula <- function(x, jitter.=FALSE, ...){
    x <- x$copula
    if (jitter.){
        x <- apply(x, 2, jitter)
    }
    pairs(x, ...)
    invisible()
}

test.copula <- function(){
  fun <- function(d) apply(d,2,function(x)(1:n)[rank(x)])/(1+n)
  n <- 200

  u2 <- cbind(sample(n),sample(n))
  d2 <- fun(u2)

  u3 <- cbind(sample(n),sample(n),sample(n))
  d3 <- fun(u3)

  checkEqualsNumeric(d2,copula(u2)$copula,msg="copula: 2 dimensional")
  checkEqualsNumeric(d3,copula(u3)$copula,msg="copula: 3 dimensional")
  checkException(copula(TRUE),msg="copula: exception")
  checkException(copula("text"),msg="copula: exception")
}
