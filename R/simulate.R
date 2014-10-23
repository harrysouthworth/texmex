# Following simulate.lm, simulate new set of responses from the whole data

simulate.evmOpt <- function(object, nsim=1, seed=NULL, param=NULL, ...){
  if (!is.null(seed)) set.seed(seed)
  if (is.null(param)){
    param <- predict(object, type="lp", unique.=FALSE)
    param <- param$link[, colnames(param$link) %in% names(param$family$param)]
    if (nsim > 1)
      param <- do.call("rbind", rep(list(param), nsim))
  }

  res <- unname(object$family$rng(nrow(param), param, object))
  if (nsim > 1)
    matrix(res, byrow=FALSE, ncol=nsim)
  else
    res
}

simulate.evmSim <- function(object, nsim=1, seed=NULL, ...){
  if (!is.null(seed)) set.seed(seed)

  res <- list()
  for (i in 1:nsim){ # Do this here to stop reuse of same parameters
    # permute parameters to avoid reusing the same ones when nsim > 1
    object$param <- object$param[sample(nrow(object$param)), ]
    
    param <- predict(object, type="lp", unique.=FALSE)
    param <- param$link[, colnames(param$link) %in% names(param$family$param)]
    
    res[[i]] <- simulate(object$map, nsim=1, param=param, seed=NULL)
  }
  do.call("cbind", res)
}

simulate.evmBoot <- function(object, nsim=1, seed=NULL, ...){
  if (!is.null(seed)) set.seed(seed)

  res <- list()
  for (i in 1:nsim){
    # permute parameters to avoid reusing the same ones when nsim > 1
    object$replicates <- object$replicates[sample(nrow(object$replicates)), ]
    
    param <- predict(object, type="lp", unique.=FALSE)
    param <- param$link[, colnames(param$link) %in% names(param$family$param)]
    
    res[[i]] <- simulate(object$map, nsim=1, param=param, seed=NULL)
  }
  do.call("cbind", res)
}
