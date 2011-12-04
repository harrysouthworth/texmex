# Author: Harry Southworth
# Date: 2011-11-25
# Purpose: Create a predict method for objects of class gpd and bgpd that
#          returns parameters, return levels or (maybe) return periods, 
#          depending on arguments given.
#
# predict.gpd
# predict.bgpd
# predict.bootgpd

################################################################################



################################################################################
## Need methods for gpd, bgpd, bootgpd.

################################################################################
## gpd

predict.gpd <-
    # Get predictions for a gpd object. These can either be the linear predictors
    # or return levels.
function(object, newdata=NULL, type=c("return level", "link"), se.fit=FALSE, ci.fit=FALSE, M=1000){
    theCall <- match.call()
    
    type <- match.arg(type)
        
    res <- switch(type,
                  "return level" = rl.gpd(object, M, newdata),
                  "link" = predict.link.gpd(object, newdata, se.fit, ci.fit)
                  )
    res <- list(rl = res, call = theCall)
    oldClass(res) <- "returnLevel"
    res
}

## Linear predictor functions for GPD

predict.link.gpd <- function(object, newdata, se.fit, ci.fit){



}

## Return level functions for GPD

## Reversing arguments M and newdata for anyone who wants to call these functions
## directly

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



