\name{hist.gpd}
\alias{hist.gpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hist.gpd(object, xlab, ylab, main)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{xlab}{ ~~Describe \code{xlab} here~~ }
  \item{ylab}{ ~~Describe \code{ylab} here~~ }
  \item{main}{ ~~Describe \code{main} here~~ }
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
function( object , xlab , ylab , main ){
	a = object$mle
	a[ 1 ] = exp( a[ 1 ] )
	u = object$threshold
	dat = object$data

	h <- hist(dat, plot = F)
	x <- seq(u, max(h$breaks), length = 100)
	y <- dgpd( x , a[ 1 ] , a[ 2 ], u )
	
	if ( missing( xlab ) || is.null( xlab ) ) xlab = "x"
	if ( missing( ylab ) || is.null( ylab ) ) ylab = ""
	if ( missing( main ) || is.null( main ) ) main = "Histogram and density"
	
	hist( dat, prob = T, ylim = c(0, max(y)),
		  xlab = xlab, ylab = ylab, 
		  main = main
		 )
	lines(x, y, col = 4)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
