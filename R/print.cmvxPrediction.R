`print.cmvxPrediction` <-
function( x, ... ){
	if ( is.R() ) stdev <- function( x ) sqrt( var( x ) )

	print( x$call, ... )
	cat( "\nResults from", length( x$replicates ), "bootstrap runs.\n" )
# XXX Need something about effective samples here

	cv <- names( x$data$simulated )[ 1 ]
	cat( paste( "\nConditioned on ", cv, " being above its ", 100*x$pqu, "th percentile.\n\n", sep = "" ) )
	
	res <- t( sapply( x$replicates , function ( x ) apply( x, 2, mean ) ) )

	m <- apply( res, 2, mean )
	s <- apply( res, 2, stdev )
	res <- matrix( c( m, s ), byrow=TRUE, nrow=2 )
	dn <- paste( "E(", dimnames( x$replicates[[ 1 ]] )[[ 2 ]] ,"|", cv , ">Q",100*x$pqu,")", sep="" )
	dimnames( res ) <- list( c( "mean", "se" ), dn )
	print( res, ... )
		
	invisible( res )
}

