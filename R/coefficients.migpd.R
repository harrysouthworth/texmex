`coefficients.migpd` <-
function( object, ... ){
	co <- sapply( object$models, function( x ) x$coefficients )
	co <- rbind( object$mth, object$mqu, co )
	dimnames( co ) <- list( c( "Threshold", "P(X < threshold)", "sigma", "xi" ) ,
							names( object$models ) )
	co[ 3, ] <- exp( co[ 3 , ] )
	co
}
