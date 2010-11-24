`hist.gpd` <-
function( x , xlab, ylab, main, ... ){
	a <- x$coefficients
	a[ 1 ] <- exp( a[ 1 ] )
	u <- x$threshold
	dat <- x$y

	h <- hist(dat, plot = FALSE)
	x <- seq(u, max(h$breaks), length = 100)
	y <- dgpd( x , a[ 1 ] , a[ 2 ], u )

	if ( missing( xlab ) || is.null( xlab ) ) xlab = "x"
	if ( missing( ylab ) || is.null( ylab ) ) ylab = ""
	if ( missing( main ) || is.null( main ) ) main = "Histogram and density"
	
	hist( dat, prob = TRUE, ylim = c(0, max(y)),
		  xlab=xlab, ylab=ylab, main=main, ...
		 )
	lines(x, y, col = 4)
	rug(dat)
	invisible()
}

