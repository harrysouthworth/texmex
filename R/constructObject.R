addCoefficients <-
    # Add named coefficients to object returned by optim
function(o){
    o$coefficients <- o$par
    o$par <- NULL
    
    nms <- unlist(lapply(names(modelData$D),
                         function(x){
                             paste(x, ": ", colnames(modelData$D[[x]]), sep = "")
                         } ) )
    names(o$coefficients) <- nms
    o
}

addCovariance <- function(o, family, cov){
    if (cov == "numeric" | is.null(family()$info)) {
      o$cov <- solve(o$hessian)
    }
    else if (cov == "observed") {
      o$cov <- solve(family()$info(o))
    }
    else {
      stop("cov must be either 'numeric' or 'observed'")
    }

    o$hessian <- NULL
    o
}

constructEVM <- function(o, family, th, rate, prior, modelParameters, call,
                         data, priorParameters){
    o$threshold <- th
    o$rate <- rate
    o$penalty <- prior
    o$coefficients <- addCoefficients(o)
    o$formulae <- modelParameters
    o$call <- call
    o$data <- data
    o$residuals <- family()$resid(o)
    o$priorParameters <- priorParameters
    o$loglik <- -o$value
    o$value <- o$counts <- NULL
    o$cov <- addCovariance(o, family, cov)
    o$se <- sqrt(diag(o$cov))

    oldClass(o) <- 'evm'
    o
}
