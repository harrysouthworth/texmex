`ppgpd` <-
function( object , nsim = 1000, alpha = .050,
		  xlab, ylab,  main, # labels and titles
		  pch=1, col = 2, cex = .75, linecol = 4 ,
		  cicol = 0, polycol=15, smooth = TRUE ){
	a <- object$coefficients
	a[ 1 ] <- exp( a[ 1 ] )
	u <- object$threshold
	dat <- object$y

	if ( missing( xlab ) || is.null( xlab ) ) { xlab <- "Empirical" }
	if ( missing( ylab ) || is.null( ylab ) ) { ylab <- "Model" }
	if ( missing( main ) || is.null( main ) ) { main <- "Probability Plot" }

	x <- ( 1:length( dat ) ) / length( dat )

	# If doing the envelope, simulate, sort and get the quantiles
	if ( nsim > 0 ){
		n <- length( dat )
		sim <- matrix( rgpd( nsim * n, a[ 1 ], a[ 2 ], u ), ncol = nsim )
		sim <- apply( sim, 2, sort )
		sim <- apply( sim, 1, pgpd, sigma = a[ 1 ], xi = a[ 2 ] , u = u )
		sim <- apply( sim, 2, quantile, prob = c( alpha/2, 1 - alpha/2 ) )
	}
    else {
        sim <- NULL
    }

	ylimits <- if ( nsim == 0 ) NULL else range( sim )

	oldpar <- par( pty = "s" ); on.exit( oldpar )
	plot( x,
		  pgpd( sort( dat ), a[ 1 ] , a[ 2 ] , u ),
		  xlab = xlab, ylab = ylab,
		  main = main, 
		  type = "n"
		 )
	# If doing the envelope, plot it before putting the data on 
	if ( nsim > 0 ){
		if ( smooth & length( x ) > 2 ) {
			splo <- spline( x, sim[ 1 , ] , 200 )
			sphi <- spline( x, sim[ 2 , ] , 200 )
			if ( polycol != 0 )
				polygon( c( splo$x, rev( sphi$x ) ),
					     c( splo$y, rev( sphi$y ) ),
					     col = polycol ,
						border = FALSE
					    )
			lines( splo$x, splo$y, col = cicol )
			lines( sphi$x, sphi$y, col = cicol )
		}
		else{
			if ( polycol != 0 )
				polygon( c( x, rev( x ) ), c( sim[ 1, ], rev( sim[ 2, ] ) ),
						col=polycol, border=FALSE )
			lines( x, sim[ 1, ], col = cicol )
			lines( x, sim[ 2, ] , col = cicol )
		}
	}


	abline( 0, 1, col = linecol )
	points( x,
		    pgpd( sort( dat ), a[ 1 ] , a[ 2 ] , u ),
			pch = pch , col = col, cex = cex
		  )
	box()
	invisible( sim )
}

