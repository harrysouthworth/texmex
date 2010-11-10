`bgpd` <-
function(y, data, th, qu, phi = ~ 1, xi = ~ 1, prior="gaussian",
         priorParameters = NULL,
		     iter = 10500, burn=500, thin = 1 ,
		     jump.const, # Scale the covar of the jumping dist
		     start, trace = 1000,
         verbose = TRUE){

  require(mvtnorm,quietly=verbose)
  
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
            stop("wrong number of location parameters in prior (doesn't match phi and xi formulas)")
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

  if ( missing( start ) ) 
    res[ 1 , ] <- coef( mod )
	else 
    res[ 1 , ] <- start

  if (!exists(".Random.seed")) runif(1) 
	seed <- .Random.seed # Retain and add to output

	for( i in 2:iter ){
    if( verbose){
      if( i %% trace == 0 ) cat(i, " steps taken\n" )
    }
		if ( exists("is.R") && is.function(is.R) && is.R() ){
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
    if(verbose) {
      cat( "Acceptance rate:", round( acc , 3 ) , "\n" )
    }
	}
	
	res <- list( call=theCall, threshold=u , map = mod ,
					burn = burn, thin = thin, 
					chains=res ,
					X.phi = X.phi, X.xi = X.xi,
					acceptance=acc, seed=seed 
					)
  oldClass( res ) <- "bgpd"
  res <- thinAndBurn(res)
	invisible( res )
}

#*************************************************************
#*************************************************************
#*************************************************************
# TESTING ROUTINE

test(bgpd) <- function(){
  require(svUnit,quiet=TRUE)
  
  postSum <- function(x){
    t(apply(x$param, 2, function(o){ c(mean=mean(o), se=sd(o)) }))
  } 
  
#************************************************************* 
# 4.1. Test reproducibility

  set.seed(20101110)
  save.seed <- .Random.seed

  set.seed(save.seed)
  bmod <- bgpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7), iter=1000,verbose=FALSE)
  
  set.seed(save.seed)
  bmod2 <- bgpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7), iter=1000,verbose=FALSE)
  
  checkEqualsNumeric(bmod$param,bmod2$param)

  set.seed(bmod$seed)
  bmod3 <- bgpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7), iter=1000,verbose=FALSE)
  checkEqualsNumeric(bmod$param, bmod3$param)

#*************************************************************  
# 4.2. Logical test of burn-in

  checkEqualsNumeric(nrow(bmod$chains) - bmod$burn, nrow(bmod$param))

  iter <- sample(500:1000,1)
  burn <- sample(50,1)
  bmod2 <- bgpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7),
                iter=iter, burn=burn,verbose=FALSE)

  checkEqualsNumeric(iter-burn, nrow(bmod2$param))

#*************************************************************
# 4.3. Logical test of thinning

  thin <- 0.5
  iter <- 1000
  bmod <- bgpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7),
               iter=iter, thin = thin,verbose=FALSE)

  checkEqualsNumeric((nrow(bmod$chains) - bmod$burn) * thin, nrow(bmod$param))

  thin <- 2
  iter <- 1000
  bmod <- bgpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7),
               iter=iter, thin = thin,verbose=FALSE)

  checkEqualsNumeric((nrow(bmod$chains) - bmod$burn) / thin, nrow(bmod$param))

#*************************************************************  
## Checks of bgpd. Centres of distributions and standard deviations 
## should be similar to those coming out of gpd with the same penalty.
## Posterior means are not the same as MAPs and SEs from Hessian are
## approximations, so allow some tolerance.

  tol <- 0.01
  tol.se <- 0.2
# 4.4. Compare MAP and posterior summaries for simple model

  mod <- gpd(ALT_M, data=liver, qu=.7)
  bmod <- bgpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7),verbose=FALSE)
  
  checkEqualsNumeric(coef(mod), postSum(bmod)[,1], tolerance=tol)
  checkEqualsNumeric(mod$se, postSum(bmod)[,2], tolerance=tol.se)

#*************************************************************
# 4.5. Covariates in phi
  tol <- 0.01
  tol.se <- 0.2

  require(MASS,quiet=TRUE) # Use MASS because is ships with R and has no dependencies
              # that don't ship with R.
  liver <- liver
  liver$ndose <- as.numeric(liver$dose)

# function rlm : Fit a linear model by robust regression using an M estimator:

  rmod <- rlm(log(ALT_M) ~ log(ALT_B) + dose, data=liver) 
  liver$resids <- resid(rmod)

  mod <- gpd(resids, data=liver, qu=.7, phi = ~ ndose)

  bmod <- bgpd(resids, data=liver, th=quantile(liver$resids, .7),
               phi = ~ ndose,verbose=FALSE)
              
  checkEqualsNumeric(coef(mod), postSum(bmod)[,1], tolerance=tol)
  checkEqualsNumeric(mod$se, postSum(bmod)[,2], tolerance=tol.se)

#*************************************************************
# 4.6. Covariates in xi
  tol <- 0.02
  tol.se <- 0.2

  mod <- gpd(resids, data=liver, qu=.7, xi = ~ ndose)
  bmod <- bgpd(resids, data=liver, th=quantile(liver$resids, .7),
               xi = ~ ndose,verbose=FALSE)

  checkEqualsNumeric(coef(mod), postSum(bmod)[,1], tolerance=tol)
  checkEqualsNumeric(mod$se, postSum(bmod)[,2], tolerance=tol.se)  

#*************************************************************
# 4.7. Covariates in xi and phi
  tol <- 0.02
  tol.se <- 0.3

  mod <- gpd(resids, data=liver, qu=.7,
             xi = ~ ndose, 
             phi = ~ ndose)
  bmod <- bgpd(resids, data=liver, th=quantile(liver$resids, .7),
               xi = ~ ndose, 
               phi = ~ ndose,verbose=FALSE)
               
  checkEqualsNumeric(coef(mod), postSum(bmod)[,1], tolerance=tol)
  checkEqualsNumeric(mod$se, postSum(bmod)[,2], tolerance=tol.se)  
 
}