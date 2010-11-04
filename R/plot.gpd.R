`plot.gpd` <-
function( object, main=rep(NULL,4), xlab=rep(NULL,4), nsim=1000, alpha=.05 ){
    if ( !missing( main ) ){
        if ( length( main ) != 1 & length( main ) != 4 ){
            stop( "main should have length 1 or 4" )
        }
        else if ( length( main ) == 1 ){ main <- rep( main, 4 ) }
    }
    
    if (ncol(object$X.phi) == 1 && ncol(object$X.xi) == 1){
#        n <- ncol(object$X.phi)
#        x <- (1:n)/(n + 1)
        ppgpd( object, main=main[1], xlab=xlab[1], nsim=nsim, alpha=alpha )
        qqgpd( object, main=main[2], xlab=xlab[2], nsim=nsim, alpha=alpha )
        rl.gpd( object, main=main[3], xlab=xlab[3] )
        hist.gpd( object, main=main[4], xlab=xlab[4] )
    }
    else { # Covariates in the model
        sigma <- exp(coef(object)[1:ncol(object$X.phi)] %*% t(object$X.phi))
        xi <- coef(object)[(ncol(object$X.phi) + 1):length(coef(object))] %*% t(object$X.xi)
        y <- xi * (object$y - object$threshold) / sigma
        y <- 1/xi * log(1 + y) # Standard exponential
       
        object$y <- y
        object$threshold <- 0
        object$coefficients <- c(0, 0) # phi not sigma, so 0 not 1
        ppgpd( object, main=main[1], xlab=xlab[1], nsim=nsim, alpha=alpha )
        qqgpd( object, main=main[2], xlab=xlab[2], nsim=nsim, alpha=alpha )
    }

    invisible()
}


