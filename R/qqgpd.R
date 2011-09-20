qqgpd <-
function( object, nsim=1000, alpha=.050,
		  xlab, ylab, main , plot = TRUE,
		  ylim = "auto", 
		  pch= 1, col =2 , cex=.75, linecol = 4 ,
		  intcol = 0, polycol = 15 ){

	a <- object$coefficients
	a[ 1 ] <- exp( a[ 1 ] )
	u <- object$threshold
	dat <- object$y

	if ( missing( xlab ) || is.null( xlab ) ) { xlab <- "Model" }
	if ( missing( ylab ) || is.null( ylab ) ) { ylab <- "Empirical" }
	if ( missing( main ) || is.null( main ) ) { main <- "Quantile Plot" }

  ModPoints <- qgpd(ppoints(dat), a[1], a[2], u)

# If doing the envelope, simulate, sort and get the quantiles
	if ( nsim > 0 ){
		n <- length( dat )
		sim <- matrix( rgpd( nsim * n, a[ 1 ], a[ 2 ], u ), ncol = nsim )
		sim <- apply( sim, 2, sort )
		# Get the simulated MSEs
		sim.mse <- apply( sim, 2, function( x, m ) mean( (x - m)^2 ), m = ModPoints )
		sim <- apply( sim, 1, quantile, prob = c( alpha/2, 1 - alpha/2 ) )
	}
	else { sim <- NULL }
	
	# Get the quantile of the observed data MSE relative to the sims
	p <- mean( sim.mse > mean( ( sort( dat ) - ModPoints )^2 ) )

	limits <- if ( nsim == 0 ) NULL else range( sim, dat )

	if ( plot ){
		oldpar <- par( pty = "s" ); on.exit( oldpar )	
		plot( ModPoints, sort(dat), 
			  xlim = limits, ylim=limits,
			  ylab = ylab, xlab = xlab, main = main,
			  type = "n" 
			 )
		# If doing the envelope, plot it before putting the data on 
		if ( nsim > 0 ){
		  if ( polycol != 0 )
					polygon( c(ModPoints,rev(ModPoints)), c( sim[ 1, ], rev( sim[ 2, ] ) ),  col=polycol ,border=NA)
			lines( ModPoints, sim[ 1, ], col = intcol ) 
			lines( ModPoints, sim[ 2, ], col = intcol )
		}
			
		# Add the diagonal reference line and the data
		abline( 0, 1, col = linecol )
		points( ModPoints, sort(dat),
			  	pch = pch, col = col, cex=cex
			   )
		box()
	} # Close if ( plot )
	res = list( data = sort( dat ), envelope = sim, Q = p )
	invisible( res )
}

