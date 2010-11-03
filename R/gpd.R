gpd <-
function(y, data, th, qu, phi= ~1, xi= ~1,
         penalty="gaussian", priorParameters=NULL,
         maxit=10000, trace=0) {

	theCall <- match.call()

    if (class(try(y, silent=TRUE)) == "try-error"){
        if (!missing(data)){
            y <- deparse(substitute(y))
        }
    }

    if (!missing(data)){
        y <- data[, y]
        if (missing(th)){
            th <- quantile(y, qu)
        }
        X.phi <- model.matrix(phi, data)
        X.xi <- model.matrix(xi, data)
        X.phi <- X.phi[y > th,]
        X.xi <- X.xi[y > th,]
    }
    else {
        if (missing(th)){
            th <- quantile(y, qu)
        }
        X.phi <- matrix(ncol=1, rep(1, length(y)))
        X.xi <- matrix(ncol=1, rep(1, length(y)))
        X.phi <- X.phi[y > th, ]
        X.xi <- X.xi[y > th, ]
    }

    rate <- mean(y > th)

    y <- y[y > th]

    if (!is.matrix(X.phi)){
        X.phi <- matrix(X.phi, ncol=1)
    }
    if (!is.matrix(X.xi)){
        X.xi <- matrix(X.xi, ncol=1)
    }

	  if ( length( y ) == 0 ){ stop( "No observations over threshold" ) }
      if ( casefold( penalty ) %in% c( "quadratic" , "gaussian" )  ){
		if ( is.null( priorParameters ) ){
		    priorParameters <- list(rep(0, ncol(X.phi) + ncol(X.xi)),
		                              diag(rep(10^4, ncol(X.phi) + ncol(X.xi)))
		                              )
        }
      }
      else if (casefold(penalty) %in% c("lasso", "l1", "laplace")){
          if (is.null(priorParameters)){
              priorParameters <- list(rep(0, ncol(X.phi) + ncol(X.xi)),
                                      diag(rep(10^(-4), ncol(X.phi) + ncol(X.xi)))
                                      )
          }
          if (length(priorParameters) != 2 | !is.list(priorParameters) ){
              stop("For Laplace prior, priorParameters should be a list
                    of length 2, the second element of which should be
                    a diagonal matrix")
          }
          if (!is.matrix(priorParameters[[2]])){
              # Lasso penalty is a bit unintuitive. Try to help out.
              priorParameters[[2]] <- diag(rep(priorParameters[[2]], ncol(X.phi) + ncol(X.xi)))      
          }
          if (!all(wh == diag(diag(wh)))){
              warning("some off-diagonal elements of the covariance are non-zero. Only
                       the diagonal is used in penalization")
          }
      }

        if (!is.null(priorParameters)){
            dimp <- ncol(X.phi) + ncol(X.xi)
            if (length(priorParameters[[1]]) != dimp){
                stop("wrong number of parameters in prior (doesn't match phi and xi formulas)")
            }
            else if (length(diag(priorParameters[[2]])) != dimp){
                stop("wrong dimension of prior covariance (doesn't match phi and xi formulas)")
            }

        }

        gpd.lik <- function(par, y, th, X.phi, X.xi,
                            penalty = "none", priorParameters=NULL,
                            maxit=maxit, trace=trace) {
            keepsc <- par[1:ncol(X.phi)]
            keepxi <- par[ - (1:ncol(X.phi))]

            sc <- colSums(par[1:ncol(X.phi)] * t(X.phi))
            xi <- colSums(par[(ncol(X.phi) + 1):(ncol(X.phi) + ncol(X.xi))] * t(X.xi))

            y <- (y - th)/exp(sc)
            y <- 1 + xi * y

            if (min(y) <= 0) 
                l <- 10^6
            else l <- sum(sc) + sum(log(y) * (1/xi + 1)) # neg log lik
            if (casefold(penalty) == "none") 
                l <- l
            else if (casefold(penalty) %in% c("quadratic", "gaussian")) {

                p <- mahalanobis(matrix(c(keepsc, keepxi), nrow = 1), 
                  center = priorParameters[[1]], cov = priorParameters[[2]])
                l <- l + p
            }
            else if (casefold(penalty) %in% c("lasso", "l1", "laplace")){
                p <- sum(abs(c(keepsc, keepxi) - priorParameters[[1]]) * diag(priorParameters[[2]]))
                l <- l + p
            }
            else stop("penalty can be 'none', 'lasso' or 'gaussian'")
            l
        }

        start <- c(log(mean(y)), rep(.00001, -1 + ncol(X.phi) + ncol(X.xi)))

        o <- optim(par = start, fn = gpd.lik,
                   y = y, X.phi = X.phi, X.xi = X.xi, th = th, 
                   penalty = penalty, 
                   control = list(maxit = maxit, trace = trace), 
                   priorParameters = priorParameters,
		hessian =  TRUE )

      if (o$convergence != 0){
          warning("Non-convergence in gpd")
      }
      o$cov <- solve(o$hessian)
	  o$hessian <- NULL
	  o$se <- sqrt( diag( o$cov ) )
	  o$threshold <- th
	  o$penalty <- penalty
	  o$mle <- o$par
      names(o$mle) <- c(paste("phi:", colnames(X.phi)),
                        paste("xi:", colnames(X.xi)))
	  o$par <- NULL
	  o$rate <- rate
	  o$call <- theCall
	  o$y <- y
	  o$X.phi <- X.phi
	  o$X.xi <- X.xi
	  if (missing(data)) { data <- NULL }
      o$data <- data
	  o$loglik <- -gpd.lik(o$mle, y, th, X.phi, X.xi)
    o$value <- NULL # since this duplicates the $loglik value
    o$counts <- NULL # not required
	  oldClass( o ) <- "gpd"
	  o
    }

