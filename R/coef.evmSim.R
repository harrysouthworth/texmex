#' @export
coef.evmSim <- function(object, ...){
  res <- apply(o$param, 2, mean)
  names(res) <- names(o$map$coefficients)
  res
}

