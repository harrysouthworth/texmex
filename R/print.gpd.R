`print.gpd` <-
function( x , digits=max(3, getOption("digits") - 3), ... ){
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

	co <- cbind( coef( x ), x$se )
	dimnames(co) <- list(names(coef(x)) , c("Value", "SE"))
	cat( "\n\nCoefficients:\n" )
    print.default(format(co, digits=digits, ...), print.gap=2, quote=FALSE)
	cat( "\n" )
	invisible()
}
show.gpd <- print.gpd
