`summary.gpd` <-
function( object , nsim = 1000 , alpha = .050, ... ){

    if (ncol(object$X.phi) == 1 && ncol(object$X.xi) == 1){
    	env <- qqgpd(object, plot = FALSE, nsim = nsim, alpha = alpha)
    }
    else {
        x <- object
        sigma <- exp(coef(x)[1:ncol(x$X.phi)] %*% t(x$X.phi))
        xi <- coef(x)[(ncol(x$X.phi) + 1):length(coef(x))] %*% t(x$X.xi)
        y <- xi * (x$y - x$threshold) / sigma
        y <- 1/xi * log(1 + y) # Standard exponential
       
        x$y <- y
        x$threshold <- 0
        x$coefficients <- c(0, 0) # phi not sigma, so 0 not 1

        env <- qqgpd(x, plot = FALSE, nsim=nsim, alpha=alpha )
    }

	res <- list( model = object, envelope = env , nsim = nsim, alpha = alpha )
	oldClass( res ) <- "summary.gpd"
	res
}

