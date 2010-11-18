`print.summary.gpd` <-
function( x, dig = 3 , ... ){
	print( x$model , ... )
	cat( paste( x$nsim, " simulated data sets compared against observed data QQ-plot.\n", sep = "" ) )
	cat( "Quantile of the observed MSE: ", round( x$envelope$Q, ... ), "\n" )
	out = sum( x$envelope$data < x$envelope$envelope[1, ] | x$envelope$data > x$envelope$envelope[2, ] )
	level = paste( round( 100 - x$alpha*100 , ...), "%", sep = "" )
	perc = round( out / length( x$envelope$data ) * 100, dig=dig, ... )
	perc = paste( perc, "%", sep = "" )
	cat( paste( out, " observations (", perc, ") outside the ", level, " simulated envelope.\n" , sep="") )
	invisible()
}

