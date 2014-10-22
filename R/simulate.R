simulate.evmOpt <- function(object, nsim=1, seed=NULL, param=NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(param)) param <- t(coef(object))
  unname(object$family$rng(nsim, param, object))
}

simulate.evmSim <- function(object, nsim=NULL, seed=NULL, ...){
  if (!is.null(seed)) set.seed(seed)

  param <- predict(object, type="lp", unique.=FALSE)
  param <- param$link[, colnames(param$link) %in% names(param$family$param)]

  simulate(object$map, nrow(param), param=param, seed=seed)
}

simulate.evmBoot <- function(object, nsim=NULL, seed=NULL, ...){
  if (!is.null(seed)) set.seed(seed)

  param <- predict(object, type="lp", unique.=FALSE)
  param <- param$link[, colnames(param$link) %in% names(param$family$param)]

  simulate(object$map, nrow(param), param=param, seed=seed)
}
