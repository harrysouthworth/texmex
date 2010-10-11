\name{ppgpd}
\alias{ppgpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ppgpd(object, nsim = 1000, alpha = 0.05, xlab, ylab, main, pch = 1, col = 2, cex = 0.75, linecol = 4, cicol = 0, polycol = 15, smooth = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{nsim}{ ~~Describe \code{nsim} here~~ }
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
function( object , nsim = 1000, alpha = .050,
		  xlab, ylab,  main, # labels and titles
		  pch=1, col = 2, cex = .75, linecol = 4 ,
		  cicol = 0, polycol=15, smooth = TRUE ){
	a = object$mle
	a[ 1 ] = exp( a[ 1 ] )
	u = object$threshold
	dat = object$data

	if ( missing( xlab ) || is.null( xlab ) ) xlab = "Empirical"
	if ( missing( ylab ) || is.null( ylab ) ) ylab = "Model"
	if ( missing( main ) || is.null( main ) ) main = "Probability Plot"

	x = ( 1:length( dat ) ) / length( dat )

	# If doing the envelope, simulate, sort and get the quantiles
	if ( nsim > 0 ){
		n = length( dat )
		sim = matrix( rgpd( nsim * n, a[ 1 ], a[ 2 ], u ), ncol = nsim )
		sim = apply( sim, 2, sort )
		sim = apply( sim, 1, pgpd, sigma = a[ 1 ], xi = a[ 2 ] , u = u )
		sim = apply( sim, 2, quantile, prob = c( alpha/2, 1 - alpha/2 ) )
	}

	ylimits = if ( nsim == 0 ) NULL else range( sim )

	oldpar = par( pty = "s" ); on.exit( oldpar )
	plot( x,
		  pgpd( sort( dat ), a[ 1 ] , a[ 2 ] , u ),
		  xlab = xlab, ylab = ylab,
		  main = main, 
		  type = "n"
		 )
	# If doing the envelope, plot it before putting the data on 
	if ( nsim > 0 ){
		if ( smooth & length( x ) > 2 ) {
			splo = spline( x, sim[ 1 , ] , 200 )
			sphi = spline( x, sim[ 2 , ] , 200 )
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
