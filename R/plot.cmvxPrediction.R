`plot.cmvxPrediction` <-
function( x, pch=c( 1, 3 ), col=c( 2, 8 ), cex=c( 1, 1 ), ask = TRUE, ... ){

    if ( is.R() ) { 
      d <- dim( x$data$simulated )[[ 2 ]] -1
      if ( prod( par( "mfrow" ) ) < d ){
        if ( ask ) {
          op <- par(ask = TRUE)
          on.exit(par(op))
        }
      }
    }

	xdat <- x$data$real[, 1 ]
	
	for( i in 2:( dim( x$data$real )[[ 2 ]] ) ){
		ydat <- x$data$real[, i ]
		xlimits <- range( xdat , x$data$simulated[ , 1 ] )
		ylimits <- range( ydat , x$data$simulated[ , i ] )
	
		plot( xdat , ydat , xlim=xlimits , ylim=ylimits,
			  xlab = names( x$data$simulated )[ 1 ],
			  ylab = names( x$data$simulated )[ i ],
			  type = "n"
			 )
		abline( v = x$data$pth, lty=4 )
		points( x$data$simulated[, 1 ] , x$data$simulated[, i ] ,
				col=col[ 2 ] , pch=pch[ 2 ], cex=cex[ 2 ] )
		points( xdat, ydat , pch=pch[ 1 ], col=col[ 1 ], cex= cex[ 1 ] )
	}
	
	invisible()
}

