#' @export
residuals.evmOpt <- function(object,...){
    object$residuals
}

#' @export
residuals.evmSim <- function(object,...){
    object$map$residuals
}

#' @export
residuals.evmBoot <- function(object,...){
    object$map$residuals
}
