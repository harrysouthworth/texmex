`bgpd` <-
function(y, data, th, qu, phi = ~ 1, xi = ~ 1, prior="gaussian",
          priorParameters = NULL,
		  iter = 10500, burn=500, thin = 1 ,
		  jump.const, # Scale the covar of the jumping dist
		  start, trace = 1000 ){

	theCall <- match.call()

    if (class(try(y, silent=TRUE)) == "try-error"){
        if (!missing(data)){
            y <- deparse(substitute(y))
        }
    }

    if (!missing(data)){
        ys <- y
        y <- data[, y]
        if (missing(th)){
            th <- quantile(y, qu)
        }
#        y <- data[, y]
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
    
    y <- y[y > th]

    if (!is.matrix(X.phi)){
        X.phi <- matrix(X.phi, ncol=1)
    }
    if (!is.matrix(X.xi)){
        X.xi <- matrix(X.xi, ncol=1)
    }

    if (missing(jump.const)){
        jump.const <- (2.4/sqrt(ncol(X.phi) + ncol(X.xi)))^2
    }

	u <- th
	
	if ( thin < 1 ) thin <- 1 / thin
	if ( thin %% 1 > 10^(-6) ) stop("thin, or its reciprocal, should be an integer" )
	if ( burn > iter ) stop( "burn-in is longer that the whole chain" )
	
	prior <- casefold(prior)
	penalty <- prior
	# Define the LOG-priors
	if ( prior == "gaussian" ){
        if ( casefold( penalty ) %in% c( "quadratic" , "gaussian" )  ){
		    if ( is.null( priorParameters ) ){
		        priorParameters <- list(rep(0, ncol(X.phi) + ncol(X.xi)),
		                                   diag(rep(10^4, ncol(X.phi) + ncol(X.xi)))
		                              )
            }
        }
		prior <- function( param, m, co ){# log Guassian
			dmvnorm(param, m, co, log=TRUE) 
		}
	}
	else stop( "only Gaussian priors implemented" )

    if (!is.null(priorParameters)){
        dimp <- ncol(X.phi) + ncol(X.xi)
        if (length(priorParameters[[1]]) != dimp){
            stop("wrong number of parameters in prior (doesn't match phi and xi formulas)")
        }
        else if (length(diag(priorParameters[[2]])) != dimp){
            stop("wrong dimension of prior covariance (doesn't match phi and xi formulas)")
        }
    }

	gpdlik <- # Positive loglikelihood
	function(param, data, X.phi, X.xi){

		phi <- colSums(param[1:ncol(X.phi)] * t(X.phi))
		xi <- colSums(param[(ncol(X.phi) + 1):(ncol(X.phi) + ncol(X.xi))] * t(X.xi))

		n <- length(data)
		data <- xi * data / exp(phi) + 1
		if( any( data <= 0 ) ){ #| ( xi <= -.50 ) ) -(10^10)
            -(10^10)
        }
		else 
        - sum(phi) - sum(log(data) * (1/xi + 1))
	}

	# MAPs to use as starting value. 

    if (missing(data)){
        data <- data.frame(y=y)
        ys <- "y"
    }

  mod <- do.call("gpd", list(ys, data, u, phi=phi, xi=xi, penalty=penalty,
    			priorParameters=priorParameters))

  # Need to check for convergence failure here. Otherwise, end up simulating
  # proposals from distribution with zero variance in 1 dimension.
  checkNA <- any( is.na( sqrt( diag( mod$cov ) ) ) )
  if ( checkNA ){
    stop( "MAP estimates have not converged. Cannot proceed. Try a different prior" )
  }

	res <- matrix( ncol=ncol(X.phi) + ncol(X.xi), nrow=iter )

    if ( missing( start ) ) res[ 1 , ] <- coef( mod )
	else res[ 1 , ] <- start

#    x <- x[data > u]
#	data <- data[data > u] - u

        if (!exists(".Random.seed")){ runif(1) }
	seed <- .Random.seed # Retain and add to output

	for( i in 2:iter ){
		if( i %% trace == 0 ) cat(i, " steps taken\n" )
		if (is.R() ){
			prop <- c(rmvnorm( 1 , mean = coef(mod), sigma = mod$cov * jump.const ))
		}
		else {
			prop <- c(rmvnorm( 1 , mean = coef(mod), cov = mod$cov * jump.const ))
		}
		bottom <- prior( res[ i - 1 ,], priorParameters[[1]], priorParameters[[2]] ) +
		            gpdlik( res[ i - 1 , ] , y-u, X.phi, X.xi)
		top <- prior( prop, priorParameters[[1]], priorParameters[[2]]) +
		            gpdlik( prop, y-u, X.phi, X.xi)
		r <- exp( top - bottom )
		res[ i , ] <- if ( runif( 1 ) < r ) prop
					else res[ i -1 , ]
	}

	acc <- length( unique( res[ , 1 ] ) ) / iter
    if (acc < .1){ warning("Acceptance rate is low.") }
	if ( trace < iter ){
		cat( "Acceptance rate:", round( acc , 3 ) , "\n" )
	}

    if (burn > 0){
  	    param <- res[ -( 1:burn ) , ] # Remove burn-in
  	}
	wh <- 1:dim( param )[[ 1 ]] %% thin == 0
	param <- param[ wh , ]

	# Simulate probabilities of being over threshold
#	if ( mean( y > u , na.rm = TRUE ) == 1 ) p <-1
#	else
#		p <- rbeta( dim( param )[[ 1 ]], .5 + sum( mod$xdata > u , na.rm=TRUE ),
#											.5 + sum( mod$xdata <= u , na.rm = TRUE)
#					)
	
	res <- list( call=theCall, threshold=u , map = mod ,
#                    p=p,
					burn = burn, thin = thin, 
					chains=res , param=param,
					X.phi = X.phi, X.xi = X.xi,
					acceptance=acc, seed=seed 
					)
	oldClass( res ) <- "bgpd"
	invisible( res )
}

