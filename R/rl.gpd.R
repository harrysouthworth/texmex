`rl.gpd` <-
function( object, alpha = .050,
		  xlab, ylab, main ,
		  pch= 1, col =2 , cex=.75, linecol = 4 ,
		  cicol = 0, polycol = 15, smooth = TRUE ){
	a = object$mle
	a[3] <- object$threshold
	a[ 1 ] = exp( a[ 1 ] )
	u = object$threshold
	la = object$rate
#	n = object$n
	n <- length( object$y )
	npy <- 1
	mat = object$cov
	dat = object$y
	# Next line is either a relic or needs updating with X.phi, X.xi
	xdat = object$xdat
	
	a <- c(la, a)
	eps <- 1e-006
	a1 <- a2 <- a3 <- a
	a1[1] <- a[1] + eps
	a2[2] <- a[2] + eps
	a3[3] <- a[3] + eps
	jj <- seq(0, 3.75 + log10(npy), by = 0.1)
	m <- unique( c(1/la, 10^jj) )
	q <- qgpd2(m, a[ 2 ], a[ 3 ], u, la)
	d1 <- (qgpd2( m, a1[ 2 ] , a1[ 3 ], u, la ) - q)/eps
	d2 <- (qgpd2( m, a2[ 2 ] , a2[ 3 ], u, la) - q)/eps
	d3 <- (qgpd2( m, a3[ 2 ] , a3[ 3 ], u, la) - q)/eps
	d <- cbind(d1, d2, d3)

	mat <- matrix(c((la * (1 - la))/n, 0, 0, 0, mat[1, 1], mat[1, 2], 0, 
		mat[2, 1], mat[2, 2]), nc = 3)
	q.form <- function (d, m) {
	    t(as.matrix(d)) %*% m %*% as.matrix(d)
	}
	v <- apply(d, 1, q.form, m = mat)

	if ( missing( xlab ) || is.null( xlab ) ) xlab = "Return period"
	if ( missing( ylab ) || is.null( ylab ) ) ylab = "Return level"
	if ( missing( main ) || is.null( main ) ) main = "Return Level Plot"

	plot(m/npy, q,
		 log = "x",
		 type = "n",
		 xlim = c(1, max(m)/npy),
		 ylim = c(u, max(xdat, q[q > u - 1] + 1.96 * sqrt(v)[q > u - 1])), 
		 xlab = xlab, ylab = ylab, main = main
		)
	# Do polygon and CI lines
	if ( smooth & length( xdat ) > 2 ) {
		splo = spline( m[q > u - 1]/npy, q[q > u - 1] - qnorm(1 - alpha/2) * sqrt(v)[q > u - 1] , 200 )
		sphi = spline( m[q > u - 1]/npy, q[q > u - 1] + qnorm(1 - alpha/2) * sqrt(v)[q > u - 1] , 200 )
		if ( polycol != 0 )
			polygon( c( splo$x, rev( sphi$x ) ),
				     c( splo$y, rev( sphi$y ) ),
				     col = polycol ,
					border=FALSE
				    )
		lines( splo$x, splo$y, col = cicol )
		lines( sphi$x, sphi$y, col = cicol )
	}
	else{
		if ( polycol != 0 )
			polygon( c( m[q > u - 1]/npy, rev( m[q > u - 1]/npy ) ),
			 		 c( q[q > u - 1] - qnorm(1 - alpha/2) * sqrt(v)[q > u - 1],
			 		 	rev( q[q > u - 1] + qnorm(1 - alpha/2) * sqrt(v)[q > u - 1] ) ),
			 		 col=polycol,
					border = FALSE
			 		)
		lines( m[q > u - 1]/npy,
			   q[q > u - 1] + qnorm(1 - alpha/2) * sqrt(v)[q > u - 1],
		   	   col = cicol
		  	 )
		lines( m[q > u - 1]/npy,
			   q[q > u - 1] - qnorm(1 - alpha/2) * sqrt(v)[q > u - 1],
		   	   col = cicol 
		     )
	}
	
	lines( m[q > u - 1]/npy, q[q > u - 1],
		   col = linecol[ 1 ] )
	
	nl <- n - length(dat) + 1
	sdat <- sort(xdat)
	points((1/(1 - (1:n)/(n + 1))/npy)[sdat > u], sdat[sdat > u], pch= pch, col = col )	
	box()
	invisible()
}

