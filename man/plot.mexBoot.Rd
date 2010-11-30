\name{plot.mexBoot}
\alias{plot.mexBoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.mexBoot(x, which = "gpd")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{which}{ ~~Describe \code{which} here~~ }
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
function( x , which = "gpd" ){


	if ( casefold( which ) == "gpd" ) which <- 1
	else which <- 2
	
	d2 <- dim( x$boot[[ 1 ]][[ which ]] )
	
	x <- x$boot
	co <- unlist( lapply( x , function( z, wh ) z[[ wh ]], wh=which ) )
	co <- array( co, dim = c( d2[ 1 ] , d2[ 2 ] , length( co ) / prod( d2 ) ) )
	
	# This is rubbish - need titles and x-axis labels
	fun <- function( z ) {
		hist( z , prob=TRUE, eraseoutline= FALSE, col=6 )
		lines( density( z, n=100 ) )
		invisible()
	}
	apply( co , 1:2, fun )
	invisible()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
