`plot.gpd` <-
function( object, main=rep(NULL,4), xlab=rep(NULL,4) ){
	if ( !missing( main ) ){
		if ( length( main ) != 1 & length( main ) != 4 ){
			stop( "main should have length 1 or 4" )
		}
		else if ( length( main ) == 1 ){ main <- rep( main, 4 ) }
	}
	
	n <- length(object$data)
	x <- (1:n)/(n + 1)
	ppgpd( object, main=main[1], xlab=xlab[1] )
	qqgpd( object, main=main[2], xlab=xlab[2] )
	rl.gpd( object, main=main[3], xlab=xlab[3] )
	hist.gpd( object, main=main[4], xlab=xlab[4] )
	invisible()
}

