`show.gpd` <-
function( x , dig = 3 ){
	cat( "Call: " )
	print( x$call )
	if ( is.null( x$penalty ) ) cat( "\nModel fit by maximum likelihood.\n" )
	else if ( is.list( x$penalty ) ){
		if ( x$penalty$type == 2 ) type = "Quadratic" 
		else type = "Absolute" 
		cat( "\nModel fit by penalized maximum likelihood.\n" )
		cat( type, " penalties on sigma and xi:\n" )
		cat( "Penalty parameters\n" )
		cat( "sigma: ", x$penalty$lambda.sigma, "\t\txi: ", x$penalty$lambda.xi )
		cat("\n" )
	}
	else {
		cat( "\nModel fit by penalized maximum likelihood.\n" )
		if ( casefold( x$penalty ) == "jeffreys" )
			cat( "Jeffreys' prior was used to penalize.\n" )
		else if( casefold( x$penalty ) == "coles-dixon" )
			cat( "Coles-Dixon prior was used to penalize.\n" )
		else if ( casefold( x$penalty ) == "gaussian" )
			cat( "A Gaussian prior was used to penalize.\n" )
	}
	if ( x$conv == 0 ) conv <- TRUE
	else conv <- FALSE
	cat( "\nConvergence:\t", conv, "\n" )
	cat( "\nThreshold:\t\t", round( x$threshold, dig), "\n" )
	cat( "Rate of excess:\t", round( x$rate, dig ), "\n" )
    cat("log-likelihood:\t", round(x$loglik, dig), "\n")
    cat("AIC:\t", round(-2 * x$loglik + 2 * length(coef(x)), dig), "\n")
	co <- cbind( coef( x ), x$se )
	dimnames(co) <- list(names(coef(x)) , c("Value", "SE"))
	cat( "\nCoefficients:\n" )
	print( round( co, dig=dig ) )
	cat( "\n" )
}

