\name{gpd.pwm}
\alias{gpd.pwm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gpd.pwm(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
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
function( x ){
	m = mean( x )
	x = sort( x )
	n = length( x )
	p = ( 1:n - .35 ) / n
	s = mean( ( 1 - p ) * x )
	xi = ( m / (m - 2 * s ) ) -2
	b = 2 * m * s / ( m - 2 * s )
	if ( x[ n ] < b / xi ){
		res = c( sigma = b, xi = xi )
		attr( res , "type" ) = "PWM"
	}
	else {
		res = c( sigma = b, xi = b/x[ n ] )
		attr( res, "type" ) = "Hybrid"
	}
	res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
