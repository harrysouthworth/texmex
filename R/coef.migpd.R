`coef.migpd` <-
function( object, ... ){
	co <- sapply( object$models, function( x ) x$par )
	co <- rbind( object$th, object$qu, co )
	dimnames( co ) <- list( c( "Threshold", "P(X < threshold)", "sigma", "xi" ) ,
							names( object$models ) )
	co[ 3, ] <- exp( co[ 3 , ] )
	co
}

