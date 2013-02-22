endPoint <- function(y,verbose=TRUE,...){
  UseMethod("endPoint", y)
}

endPoint.evmOpt <- function(y, verbose=TRUE,...){
#  fittedScale <- c(fittedGPDscale(y))
#  fittedShape <- c(fittedGPDshape(y))

  p <- texmexMakeParams(coef(y), y$data$D)
  endpoint <- y$family$endpoint

  negShape <- p[, ncol(p)] < 0

  if(any(negShape)){
    UpperEndPoint <- endpoint(p, y)
    UpperEndPoint[!negShape] <- Inf
    if(verbose){
      o <- unique(cbind(y$data$D[['xi']], p))
      print(signif(o,...))
    } else {
      invisible(unique(UpperEndPoint))
    }
  } else {
    Inf
  }
}

endPoint.evmBoot <- endPoint.evmSim <- function(y,verbose=TRUE,...){
  endPoint(y$map,verbose=verbose,...)
}

test.endPoint <- function(){
  set.seed(1)
  fit <- evm(rnorm(100),th=0.3)
  ep <- endPoint(fit,verbose=FALSE)
  co <- coef(fit)
  th <- fit$thresh
  checkEqualsNumeric(ep, th-exp(co[1])/co[2], msg="endPoint: check calc for evmOpt single covariate")

  set.seed(1)
  fit <- evm(rnorm(100),th=0.3,method="simulate",verbose=FALSE,iter=1500, thin=0)
  ep <- endPoint(fit,verbose=FALSE)
  co <- coef(fit$map)
  th <- fit$map$thresh
  checkEqualsNumeric(ep, th-exp(co[1])/co[2], msg="endPoint: check calc for evmSim single covariate")
}
