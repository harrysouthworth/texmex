AIC.gpd <- function(object, ..., k=2){
  -2*object$loglik + k*length(coef(object))

}
