`show.cmvxBoot` <-
function( x , na.action=na.exclude ){
	print( x$call )
	cat( paste( "\n", x$B, " bootstrap samples created.\n" , sep = "" ) )
		
	co <- coef( x, which = 1 )
	eff <- attributes( co )$"Effective samples"
	
#	cat( paste( eff, "effective samples.\n" ) )
	
	cat( "\nMarginal GPD bootstrap parameter estimates:\n" )
	print( co[[ 1 ]] )
	
	co <- coef( x, which = 2 )[[ 1 ]]
	cat( "\nDependence structure bootstrap parameter estimates:\n" )
	print( co )
}

