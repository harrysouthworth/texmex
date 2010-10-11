\name{print.summary.gpd}
\alias{print.summary.gpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.summary.gpd(x, dig = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{dig}{ ~~Describe \code{dig} here~~ }
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
function( x, dig = 3 ){
	print( x$model , dig = dig )
	cat( paste( x$nsim, " simulated data sets compared against observed data QQ-plot.\n", sep = "" ) )
	cat( "Quantile of the observed MSE: ", round( x$envelope$Q, dig=dig ), "\n" )
	out = sum( x$envelope$data < x$envelope$envelope[1, ] | x$envelope$data > x$envelope$envelope[2, ] )
	level = paste( round( 100 - x$alpha*100 , dig = dig), "\%", sep = "" )
	perc = round( out / length( x$envelope$data ) * 100, dig = dig )
	perc = paste( perc, "\%", sep = "" )
	cat( paste( out, " observations (", perc, ") outside the ", level, " simulated envelope.\n" , sep="") )
	invisible()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
