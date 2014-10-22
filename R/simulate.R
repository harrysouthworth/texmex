simulate.evmOpt <- function(object, nsim=1, seed=NULL, param=NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(param)) param <- t(coef(object))
  unname(object$family$rng(nsim, param, object))
}

simulate.evmSim <- function(object, nsim=NULL, seed=NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  simulate(object$map, nrow(object$param), param=object$param, seed=seed)
}

simulate.evmBoot <- function(object, nsim=NULL, seed=NULL, ...){
  if (!is.null(seed)) set.seed(seed)
  simulate(object$map, nrow(object$replicates), param=object$replicates, seed=seed)
}
