\name{rl.gpd}
\alias{rl.gpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rl.gpd(object, alpha = 0.05, xlab, ylab, main, pch = 1, col = 2, cex = 0.75, linecol = 4, cicol = 0, polycol = 15, smooth = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{alpha}{ ~~Describe \code{alpha} here~~ }
  \item{xlab}{ ~~Describe \code{xlab} here~~ }
  \item{ylab}{ ~~Describe \code{ylab} here~~ }
  \item{main}{ ~~Describe \code{main} here~~ }
  \item{pch}{ ~~Describe \code{pch} here~~ }
  \item{col}{ ~~Describe \code{col} here~~ }
  \item{cex}{ ~~Describe \code{cex} here~~ }
  \item{linecol}{ ~~Describe \code{linecol} here~~ }
  \item{cicol}{ ~~Describe \code{cicol} here~~ }
  \item{polycol}{ ~~Describe \code{polycol} here~~ }
  \item{smooth}{ ~~Describe \code{smooth} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
	n <- length( object$xdat )
	npy <- 1
	mat = object$cov
	dat = object$data
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
	    t(as.matrix(d)) \%*\% m \%*\% as.matrix(d)
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
		splo = spline( m[q > u - 1]/npy, q[q > u - 1] - 1.96 * sqrt(v)[q > u - 1] , 200 )
		sphi = spline( m[q > u - 1]/npy, q[q > u - 1] + 1.96 * sqrt(v)[q > u - 1] , 200 )
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
			 		 c( q[q > u - 1] - 1.96 * sqrt(v)[q > u - 1],
			 		 	rev( q[q > u - 1] + 1.96 * sqrt(v)[q > u - 1] ) ),
			 		 col=polycol,
					border = FALSE
			 		)
		lines( m[q > u - 1]/npy,
			   q[q > u - 1] + 1.96 * sqrt(v)[q > u - 1],
		   	   col = cicol
		  	 )
		lines( m[q > u - 1]/npy,
			   q[q > u - 1] - 1.96 * sqrt(v)[q > u - 1],
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
