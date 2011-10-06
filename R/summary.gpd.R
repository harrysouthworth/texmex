`summary.gpd` <-
function( object , nsim = 1000 , alpha = .050, ... ){

    if (ncol(object$X.phi) == 1 && ncol(object$X.xi) == 1){
    	env <- qqgpd(object, plot = FALSE, nsim = nsim, alpha = alpha)
    }
    else {
        x <- object
        x$y <- object$residuals # Standard exponential if GPD model true
        x$threshold <- 0
        x$coefficients <- c(0, 0) # phi not sigma, so 0 not 1

        env <- qqgpd(x, plot = FALSE, nsim=nsim, alpha=alpha )
    }

    co <- cbind(object$coefficients, object$se, object$coefficients / object$se)
    dimnames(co) <- list(names(coef(object)), c("Value", "SE", "t"))
    
	res <- list( model = object, coefficients=co, envelope = env , nsim = nsim, alpha = alpha )
	oldClass( res ) <- "summary.gpd"
	res
}

`print.summary.gpd` <-
    function(x, digits = 3 , ... ){
    
    co <- coef(x)
    env <- x$envelope
    nsim <- x$nsim
    alpha <- x$alpha
    
    x <- x$model

    cat( "Call: " )
    print( x$call, ... )
    if ( is.null( x$penalty ) | x$penalty=="none" ){
        cat( "\nModel fit by maximum likelihood.\n" )
    }
    else {
        cat( "\nModel fit by penalized maximum likelihood.\n" )
    }
    if ( x$conv == 0 ) conv <- TRUE
    else conv <- FALSE
    cat( "\nConvergence:\t\t")
    cat(conv)
    cat( "\nThreshold:\t\t")
    cat(format(unname(x$threshold), digits=digits, ...))
    cat( "\nRate of excess:\t\t")
    cat(format(x$rate, digits=digits, ...))
    
    cat("\n\nLog-lik.\t\tAIC\n")
    cat(format(x$loglik, digits, ...), "\t\t", format(AIC(x), digits=digits, ...))
    
    cat( "\n\nCoefficients:\n" )
    print.default(format(co, digits=digits, ...), print.gap=2, quote=FALSE)
    cat( "\n" )

    cat(paste(nsim, " simulated data sets compared against observed data QQ-plot.\n", sep = "" ))
    cat( "Quantile of the observed MSE: ", signif( env$Q, ... ), "\n")
    out <- sum(env$data < env$envelope[1, ] | env$data > env$envelope[2, ])
    level <- paste(round( 100 - alpha*100 , ...), "%", sep = "")
    perc <- round(out / length( env$data ) * 100, digits=digits, ...)
    perc <- paste(perc, "%", sep = "" )
    cat( paste( out, " observations (", perc, ") outside the ", level, " simulated envelope.\n" , sep=""))
    invisible()
}

show.summary.gpd <- print.summary.gpd
