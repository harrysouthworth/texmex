#' Compute empirical distribution function
#' 
#' Compute the empirical distribution function
#' 
#' 
#' @usage edf(x, na.last = NA)
#' @param x A numeric vector
#' @param na.last How to treat missing values. See \code{\link{rank}} for
#' details.
#' @return A vector of quantiles relating to the observations in \code{x}.
#' @author Harry Southworth
#' @seealso \code{\link{copula}}
#' @keywords univar
#' @examples
#' 
#' plot(winter$NO, edf(winter$NO))   
#' 
#' @export edf
edf <- function(x, na.last=NA){
    res <- rank(x) / (length(x) + 1)
    oldClass(res) <- "edf"
    invisible(res)
}



#' Calculate the copula of a matrix of variables
#' 
#' Returns the copula of several random variables.
#' 
#' The result is obtained by applying \code{\link{edf}} to each column of
#' \code{x} in turn.
#' 
#' Print and plot methods are available for the copula class.
#' 
#' @aliases copula plot.copula print.copula
#' @usage copula(x, na.last = NA)
#' \method{plot}{copula}(x, jitter. = FALSE, jitter.factor=1, ...)
#' \method{print}{copula}(x, ...)
#' @param x A matrix or data.frame containing numeric variables.
#' @param na.last How to treat missing values. See \code{rank} for details.
#' @param jitter. In the call to \code{plot.copula}, if \code{jitter=TRUE}, the
#' values are jittered before plotting. Defaults to \code{jitter. = FALSE}.
#' @param jitter.factor How much jittering to use. Defaults to
#' \code{jitter.factor = 1.}
#' @param ... Further arguments to be passed to plot method.
#' @return A matrix with the same dimensions as \code{x}, each column of which
#' contains the quantiles of each column of \code{x}. This object is of class
#' \code{copula}.
#' @author Harry Southworth
#' @seealso \code{\link{edf}}
#' @keywords multivariate
#' @examples
#' 
#'   D <- liver[liver$dose == "D",]
#'   Dco <- copula(D)
#'   plot(Dco)
#' 
#' @export copula
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


#' @export
print.copula <- function(x, ...){
    print(x$call)
    cat("A copula of", ncol(x$copula), "variables.\n")
    invisible(x)
}

#' @export
plot.copula <- function(x, jitter. = FALSE, jitter.factor=1, ...){
    x <- x$copula
    
    thecall <- match.call()
    jitter. <- FALSE
    if (is.element("jitter.", names(thecall))){
    	jitter. <- thecall[["jitter."]]
    }
    
	if (jitter.){
		x <- apply(x, 2, jitter, factor=jitter.factor)
	}
    pairs(x, ...)
    invisible()
}

