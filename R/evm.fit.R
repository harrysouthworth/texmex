evm.fit <- function(data, family, ...,
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

  log.lik <- family()$log.lik(data, ...)

  evm.lik <- function(par) {
    min(-log.lik(par), 1e6) + prior(par)
  }

  if (is.null(start)){
    start <- family()$start(data)
  }

   s <- unlist(lapply(data$D, function(x){ apply(x, 2, sd) }))
#   s <- c(apply(X.phi, 2, sd), apply(X.xi, 2, sd))
   s[s == 0] <- 1
#browser()
   o <- optim(par = start, fn = evm.lik,
              control = list(maxit = maxit,
              trace = trace, parscale=s),
              hessian = hessian)

    invisible(o)
}
