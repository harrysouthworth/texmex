addCoefficients <-
    # Add named coefficients to object returned by optim
function(o){
    coefficients <- o$par
    o$par <- NULL

    nms <- unlist(lapply(names(o$data$D),
                         function(x){
                             paste(x, ": ", colnames(o$data$D[[x]]), sep = "")
                         } ) )
    names(coefficients) <- nms
    coefficients
}

addCovariance <- function(o, family, cov){
    if (cov == "numeric" | is.null(family$info)) {
      cov <- solve(o$hessian)
    }
    else if (cov == "observed") {
      cov <- solve(family$info(o))
    }
    else if (cov == "sandwich") {
      if (is.null(family$sandwich)) {
          stop("sandwich estimator not implemented for this evm family")
      }
      I.inv <- solve(family$info(o))
      Sand <- family$sandwich(o)
      cov <- I.inv %*% Sand %*% I.inv
    }
    else {
      stop("cov must be either 'numeric', 'observed' or 'sandwich'")
    }

    cov
}

constructEVM <- function(o, family, ..., th, rate, prior, modelParameters, call,
                         modelData, data, priorParameters, cov){
    o$family <- family
    o$threshold <- th
    o$rate <- rate
    o$penalty <- prior
    o$data <- modelData
    o$coefficients <- addCoefficients(o)
    o$formulae <- modelParameters
    o$call <- call
    o$residuals <- family$resid(o)
    o$priorParameters <- priorParameters
    o$ploglik <- -o$value # Penalized loglik

    # Get unpenalized version
    ll <- family$log.lik(modelData, th, ...)
    o$loglik <- ll(o$coefficients)

    oldClass(o) <- 'evmOpt'

    o$cov <- try(addCovariance(o, family, cov))
    o$cov.method <- cov
    o$se <- try(sqrt(diag(o$cov)))

    o$value <- o$counts <- o$hessian <- NULL

    o$xlevels <- texmexGetXlevels(o$fo, data)

    o
}
