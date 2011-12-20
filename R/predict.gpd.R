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
                  "rl"=, "return level" = rl.gpd(object, M, newdata),
                  "link" = predict.link.gpd(object, newdata, se.fit, ci.fit, alpha)
                  )
    res
}

## Linear predictor functions for GPD

predict.link.gpd <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=.050){

    if (!is.null(newdata)){
        xi.fo <- object$call$xi
        phi.fo <- object$call$phi

        X.xi <- if (!is.null(xi.fo)){ model.matrix(as.formula(xi.fo), newdata) }
                else { matrix(1, nrow(newdata)) }
        X.phi <- if (!is.null(phi.fo)){ model.matrix(as.formula(object$call$phi), newdata) }
                 else { matrix(1, nrow(newdata)) }
    }

    else {
        X.xi <- object$X.xi
        X.phi <- object$X.phi
    }

    phi <- c(object$coefficients[1:ncol(X.phi)] %*% t(X.phi))
    xi <- c(object$coefficients[(ncol(X.phi) + 1):length(object$coefficients)] %*% t(X.xi))

    res <- cbind(phi, xi)

    if (ci.fit){
        phi.cov <- as.matrix(object$cov[1:ncol(X.phi), 1:ncol(X.phi)])
        xi.cov <- as.matrix(object$cov[(ncol(X.phi) + 1):length(object$coefficients), (ncol(X.phi) + 1):length(object$se)])
        phi.se <- sqrt(rowSums((X.phi %*% phi.cov) * object$X.phi))
        xi.se <- sqrt(rowSums((X.xi %*% xi.cov) * X.xi))

        z <- qnorm(1 - alpha/2)

        phi.lo <- phi - phi.se*z
        phi.hi <- phi + phi.se*z
        xi.lo <- xi - xi.se*z
        xi.hi <- xi + xi.se*z

        res <- cbind(res, phi.lo, phi.hi, xi.lo, xi.hi)
    } # Close if(ci.fit

    if (se.fit){
        if (!ci.fit){ # Because if ci.fit, phi.se and xi.se already exist
            phi.se <- sqrt(rowSums((X.phi %*% phi.cov) * object$X.phi))
            xi.se <- sqrt(rowSums((X.xi %*% xi.cov) * X.xi))
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

rl <- function(object, M, newdata, se.fit=FALSE, ci.fit=FALSE, alpha=.050, ...){
    UseMethod("rl")
}

gpd.delta <- function(a, m){
        # This is not exact if a prior (penalty) function is used, but
        # the CI is approximate anyway.
        
    out <- matrix(0, nrow=3, ncol=length(m))
        
    if (a[3] == 0){ # exponential case
        out[1,] <- exp(a[2]) / a[1]
        out[2,] <- exp(a[2]) * log(m * a[1])
    } else {
        out[1,] <- exp(a[2]) * m^a[3] * a[1]^(a[3] - 1)
        out[2,] <- exp(a[2]) / a[3] * ((m*a[1])^a[3] - 1) 
        out[3,] <- -exp(a[2]) / (a[3]*a[3]) * ( (m * a[1] )^a[3] - 1 ) +
                   exp(a[2]) / a[3] * (m * a[1])^a[3] * log(m * a[1])
    } 

   out
} 

rl.gpd <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=.050){
    co <- predict.link.gpd(object, newdata=newdata)
    co <- cbind(rep(object$rate, nrow(co)), co)    

    if (missing(newdata)){ co <- unique(co) }

    res <- object$threshold + (exp(co[,1]) / co[,2]) * (object$rate^co[,2] - 1)
    res <- cbind(RL=res)

    getse <- function(o, co, M){
        dxm <- t(apply(co, 1, gpd.delta, m=M))
        V <- matrix(c(object$rate * (1 - object$rate)/length(mod$y), 0, 0,
                      0, object$cov[1,],
                      0, object$cov[2,]), ncol = 3)

        # Get (4.15) of Coles, page 82, adjusted for phi = log(sigma)
        sqrt(mahalanobis(dxm, center=c(0, 0, 0), cov=V, inverted=TRUE))
    }

    if (ci.fit){
        se <- getse(object, co, M)
        lo <- res - qnorm(1 - alpha/2)*se
        hi <- res + qnorm(1 - alpha/2)*se

        res <- cbind(res, lo=lo, hi=hi)
    } # Close if (ci.fit

    if (se.fit){
        if (!ci.fit){
            se <- getse(object, co, M)
        }
        res <- cbind(res, se=se)
    }

    res
}

################################################################################
## bgpd

predict.bgpd <- function(object, newdata=NULL, type=c("return level", "link"), M=1000){
    theCall <- match.call()
    
    type <- match.arg(type)
    
    res <- switch(type,
                  "rl" =, "return level" = rl.bgpd(object, M),
                  "lp"=, "link" = cobgpd(object, newdata)
                  )
    res <- list(rl = res, call = theCall)
    oldClass(res) <- "returnLevel"
    res
}

predict.link.bgpd <- function(object, newdata, se.fit, ci.fit){
    if (!is.null(newdata)){
        xi.fo <- object$call$xi
        phi.fo <- object$call$phi

        X.xi <- if (!is.null(xi.fo)){ model.matrix(as.formula(xi.fo), newdata) }
                else { matrix(1, nrow(newdata)) }
        X.phi <- if (!is.null(phi.fo)){ model.matrix(as.formula(object$call$phi), newdata) }
                 else { matrix(1, nrow(newdata)) }
    }

    else {
        X.xi <- object$X.xi
        X.phi <- object$X.phi
    }

    res <- cbind(X.phi, X.xi) * object$param
    oldClass(res) <- "predict.link.bgpd"
    invisible(res)
}

rl.bgpd <- function(object, M){


}

################################################################################
## bootgpd

predict.bootgpd <- function(object, type=c("return level", "link"), M=1000){
    theCall <- match.call()

    type <- match.arg(type)

    res <- switch(type,
                  "rl"=, "return level" = rl.bgpd(object, M),
                  "lp"=, "link" = cobgpd(object, newdata)
                  )
    res <- list(rl = res, call = theCall)
    oldClass(res) <- "returnLevel"
    res

}
predict.link.bootgpd <- function(object, newdata, se.fit, ci.fit){
    if (!is.null(newdata)){
        xi.fo <- object$call$xi
        phi.fo <- object$call$phi

        X.xi <- if (!is.null(xi.fo)){ model.matrix(as.formula(xi.fo), newdata) }
                else { matrix(1, nrow(newdata)) }
        X.phi <- if (!is.null(phi.fo)){ model.matrix(as.formula(object$call$phi), newdata) }
                 else { matrix(1, nrow(newdata)) }
    }

    else {
        X.xi <- object$X.xi
        X.phi <- object$X.phi
    }

    res <- cbind(X.phi, X.xi) * object$replicates
    oldClass(res) <- "predict.link.bootgpd"
    invisible(res)

}

rl.bootgpd <- function(object, M){


}



