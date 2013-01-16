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
    data <- object$data$D
    cov <- v <- vector('list', length=length(data))

    # First get covariance matrix for each main parameter (e.g. phi=a'X1, xi=b'X2)
    wh <- 1
    for (i in 1:length(cov)){
        which <- wh:(wh -1 + ncol(data[[i]]))
        cov[[i]] <- as.matrix(object$cov[which, which])
        # cov[[]i] contains the block of the full covariance which relates to parameter[i]

       # Get the variance of the linear predictors
        v[[i]] <- rowSums(data[[i]] %*% cov[[i]]) * data[[i]]

        wh <- wh + ncol(data[[i]])
    }
    names(v) <- names(D)
    # Each element of v contains the variance for one linear predictor for every observation

    # We now need the off-diagonal elements of the covariance matrix. The dimensions
    # of the covariance will depend on the length of object$data$D

    getOffDiagonal <- function(k, x1, x2, cov){
        covar <- 0
        for (i in 1:ncol(x1)){
            for (j in 1:ncol(x2)){
                covar <- covar + x1[k, i] * x2[k, j] * object$cov[i, ncol(x1) + j]
            } # Close for j
        } # Close for i
        covar
    } # Close getOffDiagonal

    getCovEntry <- function(data){
        # Recursively produce off-diagonal elements of the covariance.
        # These are computed by row.
        x1 <- data[[1]]
        data[[1]] <- NULL
        n <- length(data)

        for (i in 1:n){
            res[[i]] <- sapply(1:nrow(x1), getOffDiagonal, x1, data[[i]], cov)
        }
        if (length(data) > 1){
            res <- c(res, getCovEntry(data))
        }
        else {
            res[[length(res) + 1]] <- getOffDiagonal(1:nrow(x1), x1, data[[1]], cov)
            res
        }
    } # Close getCovEntry

    co <- getCovEntry(data)

    # Now need to restructure to return a list of covariance matrices,
    # one element for every observation.
    getM <- function(i, variance, covariance){
        va <- lapply(variance, function(x){ x[[i]] })
        co <- lapply(covariance, function(x){ x[[i]] })

        res <- diag(unlist(va))
        res[upper.tri(res)] <-  res[lower.tri(res)] <- unlist(co)
        res
    }

    lapply(1:nrow(data[[1]]), getM, v, co)
}


texmexMakeCI <-
    # Compute CIs from point estimates and standard errors
function(params, ses, alpha){
    z <- qnorm(1 - alpha/2)

    lohi <- t(sapply(1:length(params), function(x, s, z){
                                           x <- x[[i]]
                                           s <- s[[i]]
                                           c(x - s*z, x + s*z)
                                       }, z=z, s=ses))
    lohi
}
