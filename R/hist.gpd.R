`hist.gpd` <-
function( x , xlab, ylab, main, ... ){
	a <- x$coefficients
	a[ 1 ] <- exp( a[ 1 ] )
	u <- x$threshold
  if(a[2] < 0){
    UpperEndPoint <- u-a[1]/a[2]
  } else {
    UpperEndPoint <- Inf
  }
	dat <- x$y

	h <- hist(dat, plot = FALSE)
	x <- seq(u, min(UpperEndPoint, max(h$breaks)), length = 100)
	y <- dgpd( x , a[ 1 ] , a[ 2 ], u )

	if ( missing( xlab ) || is.null( xlab ) ) xlab = "Data"
	if ( missing( ylab ) || is.null( ylab ) ) ylab = ""
	if ( missing( main ) || is.null( main ) ) main = "Histogram and density"
	
  breaks <- seq(from=min(dat),to=max(dat),len=nclass.Sturges(dat)+1)
	hist( dat, prob = TRUE, ylim = c(0, max(y)),
		  xlab=xlab, ylab=ylab, main=main, breaks = breaks, ...)
	lines(x, y, col = 4)
	rug(dat)
	invisible()
}

