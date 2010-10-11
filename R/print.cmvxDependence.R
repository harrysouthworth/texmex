`print.cmvxDependence` <-
function( x, ... ){
	print( x$call, ... )
	names( x$gqu ) <- dimnames( x$parameters )[[ 2 ]]
	cat( "\nThresholding quantiles for Gumbel-transformed data:\n" )
	print( x$gqu, ... )
	cat( "\nDependence structure parameter estimates:\n" )
	print( x$parameters, ... )
}

