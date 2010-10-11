\name{u2gpd}
\alias{u2gpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
u2gpd(u, p = 1, th = 0, sigma = 1, xi = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{ ~~Describe \code{u} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
  \item{th}{ ~~Describe \code{th} here~~ }
  \item{sigma}{ ~~Describe \code{sigma} here~~ }
  \item{xi}{ ~~Describe \code{xi} here~~ }
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
function( u , p = 1, th=0, sigma=1, xi= 0 )
	( ( ( 1 - u )/p )^(-xi) -1 ) * sigma/xi + th
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
