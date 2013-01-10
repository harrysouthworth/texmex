# Author: Harry Southworth
# Date: 2013-1-9
# Purpose: Split out worker functions to perform boring tasks for extreme
#          value modelling functions.
#
###########################################################################

texmexMethod <-
    # Take character string passed by user and coerce to standard format
function(method){
    method <- casefold(method)
    if (method %in% c("o", "opt", "optim", "optimize", "optimise")){
        method <- "o"
    }
    else if (method %in% c("s", "sim", "simulate")){
        method <- "s"
    }
    else {
        stop("method should be either 'optimize' or 'simulate'")
    }
    method
}

texmexPrior <-
    # Take character string(s) passed by user and coerce to standard format
function(prior, penalty, method){
    prior <- casefold(prior)
    penalty <- casefold(penalty)
    if (prior != penalty){ # User provided one or both of prior or penalty
      if (prior == 'gaussian'){ # User provided penalty
          prior <- penalty
      }
      else if (prior != 'gaussian' & penalty != 'guassian'){ # User provided both
        stop('Provide neither or one of prior and penatly, not both.')
      }
      # else prior was provided and penalty wasn't, which is ok
    }

    if (method == 's' & !is.element(prior, c('gaussian', 'cauchy'))){
      stop('Only gaussian or cauchy prior can be used when simulating from posterior.')
    }
    prior
}

texmexTrace <-
    # Get tracing frequency for optimizer and for Markov chain
function(trace, method){
    if (method == "o"){
      if (!is.null(trace)){ # trace provided by user
          otrace <- trace
      }
      else {
          otrace <- 0
      }
    } # Close if (method == "o"
    else{ # method == "s"
      otrace <- 0
      if (is.null(trace)){
         trace <- 10000
      }
    }
    c(otrace, trace)
}

texmexPrepareData <-
    # Get design matrices
function(y, data, params){
    D <- vector('list', length=length(params))
    if (!is.null(data)){
        y <- formula(paste(y, "~ 1"))
        y <- model.response(model.frame(y, data=data))

        for (i in 1:length(params)){
          D[[i]] <- model.matrix(params[[i]], data)
        }
    } # Close if(!is.null(data
    else {                                        # XXX UNTESTED CODEBLOCK XXX <---------- XXX
        for (i in 1:length(params)){
            if (length(as.character(phi)) == 2 & as.character(phi)[2] == "1"){
                D[[i]] <- matrix(ncol = 1, rep(1, length(y)))
            }
            else {
                D[[i]] <- model.matrix(params[[i]])
            }
        } # Close for
    } # Close else

    # Matrices with one column get coerced to vectors. Revert.
    D <- texmexReverseUnaskedCoercion(D)
    list(y=y, D=D)
}

texmexReverseUnaskedCoercion <-
    # R forces single column data into a vector
function(x){
    lapply(x, function(z){
                  if (!is.matrix(z)){ z <- matrix(z, ncol=1) }
                  z })
}

texmexThresholdData <- function(threshold, data){
    # Need to subset design matrices on y > th, so do those
    # first, then threshold y

    for (i in 1:length(data$D)){
        data$D[[i]] <- data$D[[i]][data$y > threshold, ]
    }

    data$D <- texmexReverseUnaskedCoercion(data$D)

    data$y <- data$y[data$y > threshold]
    if (length(data$y) == 0){
      stop("No observations above the threshold.")
    }

    data
}

texmexPriorParameters <-
    # Pre-process prior distribution parameters
function(prior, priorParameters, data){

    # Get total number of parameters
    nc <- sum(sapply(data$D, ncol))

    if (prior %in% c("quadratic", "gaussian")) {
        if (is.null(priorParameters)) {
            priorParameters <- list(rep(0, nc), diag(rep(10^4, nc)))
        }
        if (length(priorParameters) != 2 | !is.list(priorParameters)) {
            stop("For Gaussian prior or quadratic penalty, priorParameters should be a list of length 2, the second element of which should be a symmetric (covariance) matrix")
        }
    }
    else if (prior %in% c("lasso", "l1", "laplace")) {
        if (is.null(priorParameters)) {
            priorParameters <- list(rep(0, nc), diag(rep(10^(-4), nc)))
        }
        if (length(priorParameters) != 2 | !is.list(priorParameters)) {
            stop("For Laplace prior or L1 or Lasso penalty, priorParameters should be a list of length 2, the second element of which should be a diagonal (precision) matrix")
        }
        if (!is.matrix(priorParameters[[2]])) {
            priorParameters[[2]] <- diag(rep(priorParameters[[2]], nc))
        }
        if (!all(priorParameters[[2]] == diag(diag(priorParameters[[2]])))) {
            warning("some off-diagonal elements of the covariance are non-zero. Only the diagonal is used in penalization")
        }
    }

    #### If priorParameters given but of wrong dimension, kill
    if (!is.null(priorParameters)) {
        if (length(priorParameters[[1]]) != nc) {
            stop("wrong number of parameters in prior (doesn't match phi and xi formulas)")
        }
        else if (length(diag(priorParameters[[2]])) != nc) {
            stop("wrong dimension of prior covariance (doesn't match phi and xi formulas)")
        }
    }
    priorParameters
}

findFormulae <-
    # Find formulae in a call
function(call){
    wh <- sapply(call, function(x){ try(class(eval(x)), silent=TRUE) })
    wh <- names(wh)[wh == 'formula']
    as.list(call[wh])
}




