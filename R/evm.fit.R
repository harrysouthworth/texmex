evmFit <- function(data, family, ...,
                   prior="none", start=NULL,
                   priorParameters = NULL,
                   maxit = 10000, trace = 0, hessian=TRUE) {

  penFactory <- switch(prior,
                       laplace=texmex:::.make.lasso.penalty,
                       lasso=texmex:::.make.lasso.penalty,
                       l1=texmex:::.make.lasso.penalty,
                       quadratic=texmex:::.make.quadratic.penalty,
                       gaussian=texmex:::.make.quadratic.penalty,
                       none=texmex:::.make.dummy.penalty,
                       function() {stop("Bad penalty ref.")})

  prior <- penFactory(priorParameters)

  log.lik <- family$log.lik(data, ...)

  evm.lik <- function(par) {
    min(-log.lik(par), 1e6) + prior(par)
  }

  if (is.null(start)){
    if (is.null(family$start)){ start <- runif(length(s), -.1, .1) }
    else { start <- family$start(data) }
  }

   o <- optim(par = start, fn = evm.lik,
              control = list(maxit = maxit, trace = trace),
              hessian = hessian)

    invisible(o)
}
