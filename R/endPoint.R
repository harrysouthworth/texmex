endPoint <- function(y,verbose=TRUE,...){
  UseMethod("endPoint",y)
}

endPoint.gpd <- function(y,verbose=TRUE,...){
  fittedScale <- c(fittedGPDscale(y))
  fittedShape <- c(fittedGPDshape(y))
  negShape <- fittedShape < 0
  if(any(negShape)){
    UpperEndPoint <- c(y$threshold - fittedScale/fittedShape)
    UpperEndPoint[!negShape] <- Inf 
    o <- unique(cbind(y$X.xi,FittedShape=fittedShape,UpperEndPoint=UpperEndPoint))
    if(verbose){
      print(signif(o,...))
    } else {
      invisible(unique(UpperEndPoint))
    }
  } else {
    Inf
  }
}

endPoint.bgpd <- function(y,verbose=TRUE,...){
  endPoint(y$map,verbose=verbose,...)
}

fittedGPDscale <- function(o){
  exp(o$coefficients[1:ncol(o$data$D$phi)] %*% t(o$data$D$phi))
}

fittedGPDshape <- function(o){
  o$coefficients[(ncol(o$data$D$phi) + 1):length(o$coefficients)] %*% t(o$data$D$xi)
}

test.endPoint <- function(){
  set.seed(1)
  fit <- gpd(rnorm(100),th=0.3)
  ep <- endPoint(fit,verbose=FALSE)
  co <- coef(fit)
  th <- fit$thresh
  checkEqualsNumeric(ep, th-exp(co[1])/co[2], msg="endPoint: check calc for gpd single covariate")

  set.seed(1)  
  fit <- gpd(rnorm(100),th=0.3,method="simulate",verbose=FALSE,iter=1500)
  ep <- endPoint(fit,verbose=FALSE)
  co <- coef(fit$map)
  th <- fit$thresh
  checkEqualsNumeric(ep, th-exp(co[1])/co[2], msg="endPoint: check calc for bgpd single covariate")
}
