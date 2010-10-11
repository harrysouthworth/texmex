`coefficients.migpd` <-
function( x ){
	co <- sapply( x$models, function( x ) x$par )
	co <- rbind( x$th, x$qu, co )
	dimnames( co ) <- list( c( "Threshold", "P(X < threshold)", "sigma", "xi" ) ,
							names( x$models ) )
	co[ 3, ] <- exp( co[ 3 , ] )
	co
}

