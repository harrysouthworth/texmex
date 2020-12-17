#' @export
coef.evmSim <- function(object, ...){
  res <- apply(object$param, 2, mean)
  names(res) <- names(object$map$coefficients)
  res
}

