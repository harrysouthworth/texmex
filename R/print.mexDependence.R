`print.mexDependence` <-
function( x, ... ){
	print( x$call, ... )
  cat("\nConditioning on ",x$conditioningVariable," variable.\n", sep="")
	names( x$dqu ) <- dimnames( x$parameters )[[ 2 ]]
	cat( "\nThresholding quantiles for Gumbel-transformed data:\n" )
	print( x$dqu, ... )
	cat( "\nDependence structure parameter estimates:\n" )
	print( x$coefficients, ... )
}

