`summary.migpd` <-
function( object, verbose = TRUE, ... ){
	d <- dim( object$data )[ 2 ]
	conv <- sapply( object$models, function( z ) z$convergence )

	if ( sum( conv ) != 0 ) 
		conv <- ( 1:d )[ conv != 0 ]	
	else conv <- NULL
  
  co <- coefficients.migpd( object )
  
  if(verbose){
    print( object$call, ... )
    cat( "\nA collection of", d, "generalized Pareto models.\n" )
    
    if ( is.null( conv ) ) cat( "All models converged.\n" )
    else cat( "The following model(s) did not converge:", paste( conv, collapse="," ), "\n" )
    
    cat( "Penalty to the likelihood:", object$penalty )
    cat( "\n\nSummary of models:\n" )
    print( co, ... )
    cat( "\n" )
	
    if ( is.null( object$gumbel ) )
      cat( "Gumbel transformations have not yet been performed.\n" )
  }
  
	invisible( co )
}

