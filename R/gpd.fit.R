.make.quadratic.penalty <- function(priorParameters) {
  centre <- priorParameters[[1]]
  cov <- priorParameters[[2]]

  factors <- svd(as.matrix(cov))
  if (min(factors$d) == 0) {
    stop("Singular covariance matrix: quadratic penalty impossible")
  }

  stacked <- rbind(t(factors$u), t(factors$v))
  n.prior <- length(centre)
  result <- function(param) {
    delta <- param - centre
    prod <- as.numeric(stacked %*% delta)
    first <- prod[1:n.prior]
    last <- prod[(1 + n.prior):(2*n.prior)]
    as.numeric(first %*% (last / factors$d))
  }
  result
}

.make.lasso.penalty <- function(priorParameters) {
  centre <- priorParameters[[1]]
  cov <- diag(priorParameters[[2]])
  result <- function(param) {
    sum(abs(param - centre) * cov)
  }
  result
}

.make.dummy.penalty <- function(priorParameters) {
  result <- function(param) {
    0
  }
  result
}

gpdFit <- function(y, th, X.phi, X.xi, penalty="none", start=NULL,
        priorParameters = NULL, maxit = 10000, trace = 0, hessian=TRUE) {

   gpd.lik <- function(par, y, th, X.phi, X.xi, penalty) {
     n.phi <- ncol(X.phi)
     phi <- X.phi %*% par[1:n.phi]
     xi <- X.xi %*% par[(1 + n.phi):length(par)]
     tmp <- dgpd(y, exp(phi), xi, u=th, log.d=TRUE) # get the log.likelihood
     l <- -sum(tmp)
     l <- min(l, 1e6)
     l + penalty(par)
    }

   penFactory <- switch(penalty,
                        laplace=.make.lasso.penalty,
                        lasso=.make.lasso.penalty,
                        l1=.make.lasso.penalty,
                        quadratic=.make.quadratic.penalty,
                        gaussian=.make.quadratic.penalty,
                        none=.make.dummy.penalty,
                        function() {stop("Bad penalty ref.")})
   penalty <- penFactory(priorParameters)
   if (is.null(start)){
     start <- c(log(mean(y)), rep(1e-05, -1 + ncol(X.phi) + ncol(X.xi)))
   }

   s <- c(apply(X.phi, 2, sd), apply(X.xi, 2, sd))
   s[s == 0] <- 1
   o <- optim(par = start, fn = gpd.lik, y = y, X.phi = X.phi, 
              X.xi = X.xi, th = th, penalty = penalty,
              control = list(maxit = maxit,
                trace = trace, parscale=s),
              hessian = hessian)

    invisible(o)
}
