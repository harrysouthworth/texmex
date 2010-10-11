`print.migpd` <-
function( x , ... ){
	print( x$call, ... )
	cat( "\nA collection of ", dim( x$data )[ 2 ] , "generalized Pareto models.\n" )
	conv <- sapply( x$models , function( x ) x$convergence )
	if ( sum( conv ) == 0 ) 
		cat( "\nAll models converged.\n" )
	else cat( "\nWARNING: Not all models converged.\n" )

	if ( is.null( x$gumbel ) )
		cat( "Gumbel transformations have not yet been performed.\n" )

	invisible()
}

