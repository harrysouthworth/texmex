#' @export
coef.evmSim <- function(object, ...){
  object$coef
}

#' Get simulated parameters from an evm object.
#' @param object An object of class 'evmSim' or 'evmBoot'.
#' @export
param <- function(object){
  if (inherits(object, "evmSim")){
    object$param
  } else if (inherits(object, "evmBoot")){
    object$replicates
  }
}
