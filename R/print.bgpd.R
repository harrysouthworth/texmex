`print.bgpd` <-
function( x , print.seed=FALSE, ...){

	print( x$call )
	if ( print.seed ){
		cat( "Random seed:", x$seed, "\n" )
	}
	cat( "Acceptance rate: " )
	print(x$acc, ... )
	cat( "\n\nMAP estimates using", x$map$penalty, "prior:\n" )
	print( coef( x$map ) , ... )
	cat( "\nPosterior means:\n" )
	m <- c( apply( x$param, 2, mean ) )
	names( m ) = names( coef( x$map ) )
	print( m , ... )
	
	invisible()
}

