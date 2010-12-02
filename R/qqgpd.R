`qqgpd` <-
function( object, nsim=1000, alpha=.050,
		  xlab, ylab, main , plot = TRUE,
		  ylim = "auto", 
		  pch= 1, col =2 , cex=.75, linecol = 4 ,
		  cicol = 0, polycol = 15, smooth = TRUE ){

	a <- object$coefficients
	a[ 1 ] <- exp( a[ 1 ] )
	u <- object$threshold
	dat <- object$y

	if ( missing( xlab ) || is.null( xlab ) ) { xlab <- "Empirical" }
	if ( missing( ylab ) || is.null( ylab ) ) { ylab <- "Model" }
	if ( missing( main ) || is.null( main ) ) { main <- "Quantile Plot" }

	x <- qgpd( 1 - (1:length(dat)/(length(dat) + 1)), a[ 1 ], a[ 2 ] , u )

    x <- qgpd(ppoints(dat), a[1], a[2], u)

# If doing the envelope, simulate, sort and get the quantiles
	if ( nsim > 0 ){
		n <- length( dat )
		sim <- matrix( rgpd( nsim * n, a[ 1 ], a[ 2 ], u ), ncol = nsim )
		sim <- apply( sim, 2, sort )
		# Get the simulated MSEs
		sim.mse <- apply( sim, 2, function( x, m ) mean( (x - m)^2 ), m = x )
		sim <- apply( sim, 1, quantile, prob = c( alpha/2, 1 - alpha/2 ) )
	}
	else { sim <- NULL }
	
	# Get the quantile of the observed data MSE relative to the sims
	p <- mean( sim.mse > mean( ( sort( dat ) - x )^2 ) )

	ylimits <- if ( nsim == 0 ) NULL else range( sim )

	if ( plot ){
		oldpar <- par( pty = "s" ); on.exit( oldpar )	
		plot( x,
			  sort(dat),
			  ylim = ylimits,
			  ylab = ylab, xlab = xlab, main = main,
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
					polygon( c( x, rev( x ) ), c( sim[ 1, ], rev( sim[ 2, ] ) ), col=polycol )
				lines( x, sim[ 1, ], col = cicol )
				lines( x, sim[ 2, ] , col = cicol )
			}
		}
	
		# Add the diagonal reference line and the data
		abline( 0, 1, col = linecol )
		points( x,
			  	sort(dat),
			  	pch = pch, col = col, cex=cex
			   )
		box()
	} # Close if ( plot )
	res = list( data = sort( dat ), envelope = sim, Q = p )
	invisible( res )
}

