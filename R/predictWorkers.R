addCov <- function(res, X){ # used in linearPredictors.* to add covariates to columns reported in output
  if(!is.null(dim(X))){
    if(dim(X)[2] > 1){
       cov <- X[,colnames(X) != "(Intercept)"]
       res <- cbind(res,cov)
       if(is.vector(cov)) colnames(res)[dim(res)[2]] <- colnames(X)[colnames(X) != "(Intercept)"]
    }
    else {
      if( any(X != 1) ){
        res <- cbind(res, X)
      }
    }
  }
  res
}

texmexMakeParams <-
    # Take parameter vector and list of datasets to compute
    # parameters for each row of the data.
function(co, data){

    np <- length(data)
    p <- vector('list', length=np)

    wh <- 1
    for (i in 1:np){
        which <- wh:(wh -1 + ncol(data[[i]]))
        p[[i]] <- c(co[which] %*% t(data[[i]]))
        wh <- wh + ncol(data[[i]])
    }

    do.call('cbind', p)
}


texmexMakeCovariance <-
    # Get covariance matrix for each parameter and get coefficents for each row of the data
function(object){
    cov <- se <- vector('list', length=length(object$data$D))

    wh <- 1
    for (i in 1:length(res)){
        which <- wh:(wh -1 + ncol(data[[i]]))
        cov[[i]] <- as.matrix(object$cov[which, which])
        se[[i]] <- sqrt(rowSums(object$data$D[[i]] %*% cov[[i]]) * object$data$D)
    }
    list(cov=cov, se=se)
}

texmexMakeCI <-
    # Compute CIs from point estimates and standard errors
function(params, ses, alpha){
    z <- qnorm(1 - alpha/2)

    lohi <- t(sapply(1:length(params), function(x, s, z){
                                           x <- x[[i]]
                                           x <- s[[i]]
                                           c(x - s*z, x + s*z)
                                       }, z=z, s=ses))
    lohi
}
