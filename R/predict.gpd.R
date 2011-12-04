# Author: Harry Southworth
# Date: 2011-11-25
# Purpose: Create a predict method for objects of class gpd and bgpd that
#          returns parameters, return levels or (maybe) return periods, 
#          depending on arguments given.
#
################################################################################



################################################################################
## First, set up the predict wrappers

predict.gpd <- function(object, newdata, type=c("return level", "parameters"), M=1000){
    theCall <- match.call()
    
    type <- match.arg(type)
    
    if missing(newdata){
        
    
    }
    
    
    res <- switch(type,
                  "return level" = rl.gpd(object, M, newdata),
                  "parameters" = cogpd(object, newdata)
                  )
    res <- list(rl = res, call = theCall)
    oldClass(res) <- "returnLevel"
    res
}

predict.bgpd <- function(object, type=c("return level", "parameters"), M=1000){
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

################################################################################
## Return level functions for GPD

## Reversing arguments M and newdata for anyone who wants to call these functions
## directly

rl <- function(object, M, newdata, ...){
    UseMethod("rl")
}


rl.gpd <- function(object, M, newdata){
    


}

rl.bgpd <- function(object, M){


}

################################################################################
## Get parameters for GPD


