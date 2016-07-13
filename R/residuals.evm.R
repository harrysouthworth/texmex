#' @export
residuals.evmOpt <- function(object){
    object$residuals
}

#' @export
residuals.evmSim <- residuals.evmBoot <- function(object){
    object$map$residuals
}
