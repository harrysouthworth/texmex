#' @export
coef.evmSim <- function(object, ...){
  object$family$coef(object)
}
