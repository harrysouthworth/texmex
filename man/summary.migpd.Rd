\name{summary.migpd}
\alias{summary.migpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.migpd(x)
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
	d <- dim( x$data )[ 2 ]
	conv <- sapply( x$models, function( z ) z$convergence )

	if ( sum( conv ) != 0 ) 
		conv <- ( 1:d )[ conv != 0 ]	
	else conv <- NULL

	print( x$call )
	cat( "\nA collection of", d, "generalized Pareto models.\n" )
	if ( is.null( conv ) ) cat( "All models converged.\n" )
	else cat( "The following model(s) did not converge:", paste( conv, collapse="," ), "\n" )
	
	cat( "Penalty to the likelihood:", x$penalty )
	
	cat( "\n\nSummary of models:\n" )
	co <- coefficients.migpd( x ) 
	print( co )
	cat( "\n" )
	
	if ( is.null( x$gumbel ) )
		cat( "Gumbel transformations have not yet been performed.\n" )
	invisible( co )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
