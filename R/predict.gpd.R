# Author: Harry Southworth
# Date: 2011-11-25
## Purpose: Create a predict method for objects of class gpd and bgpd that
##          returns parameters, return levels or (maybe) return periods, 
##          depending on arguments given.
#
# predict.gpd
# predict.bgpd
# predict.bootgpd
# predict.link.gpd
# rl
# rl.gpd
# predict.link.bgpd
# rl.bgpd
# predict.link.bootgpd
# rl.bootgpd


################################################################################
## gpd

predict.gpd <-
    # Get predictions for a gpd object. These can either be the linear predictors
    # or return levels.
function(object, newdata=NULL, type=c("return level", "link"), se.fit=FALSE, ci.fit=FALSE, M=1000, alpha=.050){
    theCall <- match.call()
    
    type <- match.arg(type)
    
    res <- switch(type,
                  "return level" = rl.gpd(object, M, newdata),
                  "link" = predict.link.gpd(object, newdata, se.fit, ci.fit, alpha)
                  )
 #   attributes(res, "call") <- theCall
    
    res
}

## Linear predictor functions for GPD

predict.link.gpd <- function(object, newdata, se.fit, ci.fit, alpha){

    if (!is.null(newdata)){
        xi.fo <- object$call$xi
        phi.fo <- object$call$phi
#browser()
        X.xi <- if (!is.null(xi.fo)){ model.matrix(as.formula(xi.fo), newdata) }
                else { matrix(1, nrow(newdata)) }
        X.phi <- if (!is.null(phi.fo)){ model.matrix(as.formula(object$call$phi), newdata) }
                 else { matrix(1, nrow(newdata)) }
    }

    else {
        X.xi <- object$X.xi
        X.phi <- object$X.phi
   }

    phi <- c(object$coefficients[1:ncol(object$X.phi)] %*% t(object$X.phi))
    xi <- c(object$coefficients[(ncol(object$X.phi) + 1):length(object$coefficients)] %*% t(object$X.xi))

    res <- cbind(phi, xi)

    if (ci.fit){
        phi.cov <- as.matrix(object$cov[1:ncol(object$X.phi), 1:ncol(object$X.phi)])
        xi.cov <- as.matrix(object$cov[(ncol(object$X.phi) + 1):length(object$coefficients), (ncol(object$X.phi) + 1):length(object$se)])
        phi.se <- sqrt(rowSums((object$X.phi %*% phi.cov) * object$X.phi))
        xi.se <- sqrt(rowSums((object$X.xi %*% xi.cov) * object$X.xi))

        z <- qnorm(1 - alpha/2)

        phi.lo <- phi - phi.se*z
        phi.hi <- phi + phi.se*z
        xi.lo <- xi - xi.se*z
        xi.hi <- xi + xi.se*z

        res <- cbind(res, phi.lo, phi.hi, xi.lo, xi.hi)
    } # Close if(ci.fit

    if (se.fit){
        if (!ci.fit){ # Because if ci.fit, phi.se and xi.se already exist
            phi.se <- sqrt(rowSums((phi %*% phi.cov) * phi))
            xi.se <- sqrt(rowSums((xi %*% xi.cov) * xi))
        }
        res <- cbind(res, phi.se, xi.se)
    } # Close if(se.fit

    res
}

## Return level functions for GPD

## Reversing arguments M and newdata for anyone who wants to call these functions
## directly

## Will want to get return levels when using GEV rather than GPD, so make
## rl generic

rl <- function(object, M, newdata, ...){

    UseMethod("rl")
}

rl.gpd <- function(object, M, newdata){
    cat("meow\n")
    invisible()
}



################################################################################
## bgpd

predict.bgpd <- function(object, newdata=NULL, type=c("return level", "link"), M=1000){
    theCall <- match.call()
    
    type <- match.arg(type)
    
    res <- switch(type,
                  "return level" = rl.bgpd(object, M),
                  "parameters" = cobgpd(object, newdata)
                  )
    res <- list(rl = res, call = theCall)
    oldClass(res) <- "returnLevel"
    res
}

predict.link.bgpd <- function(object, newdata, se.fit, ci.fit){


}

rl.bgpd <- function(object, M){


}

################################################################################
## bootgpd

predict.bootgpd <- function(object, type=c("return level", "link"), M=1000){
    theCall <- match.call()

    type <- match.arg(type)

    res <- switch(type,
                  "return level" = rl.bgpd(object, M),
                  "parameters" = cobgpd(object, newdata)
                  )
    res <- list(rl = res, call = theCall)
    oldClass(res) <- "returnLevel"
    res

}
predict.link.bootgpd <- function(object, newdata, se.fit, ci.fit){


}

rl.bootgpd <- function(object, M){


}



