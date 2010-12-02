`print.bgpd` <-
function( x , print.seed=FALSE, ...){

	print( x$call )
	if ( print.seed ){
		cat( "Random seed:", x$seed, "\n" )
	}
	cat( "Acceptance rate: ", format(x$acc) )

	cat( "\n\nMAP estimates:\n" )
    co <- coef(x)
    map <- x$map$par
    names(map) <- names(co)

	print( map , ... )
	cat( "\nPosterior means:\n" )
	m <- coef(x)
	print( m , ... )
	
	invisible()
}

