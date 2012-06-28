

gpdFit <- function(y, th, X.phi, X.xi, penalty="none", start=NULL,
                   priorParameters = NULL,
                   maxit = 10000, trace = 0, hessian=TRUE) {

  log.lik <- .make.gpd.loglikelihood(y, th, X.phi, X.xi)

  penFactory <- switch(penalty,
                       laplace=.make.lasso.penalty,
                       lasso=.make.lasso.penalty,
                       l1=.make.lasso.penalty,
                       quadratic=.make.quadratic.penalty,
                       gaussian=.make.quadratic.penalty,
                       none=.make.dummy.penalty,
                       function() {stop("Bad penalty ref.")})
  penalty <- penFactory(priorParameters)

  gpd.lik <- function(par) {
    min(-log.lik(par), 1e6) + penalty(par)
  }

  if (is.null(start)){
    start <- c(log(mean(y)), rep(1e-05, -1 + ncol(X.phi) + ncol(X.xi)))
  }

   s <- c(apply(X.phi, 2, sd), apply(X.xi, 2, sd))
   s[s == 0] <- 1
   o <- optim(par = start, fn = gpd.lik,
              control = list(maxit = maxit,
                trace = trace, parscale=s),
              hessian = hessian)

    invisible(o)
}
