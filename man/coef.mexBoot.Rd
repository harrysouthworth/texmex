\name{coef.bootmex}
\alias{coef.bootmex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract bootstrap parameter estimates from an object of class 'bootmex' }
\description{
Extract bootstrap parameter estimates from an object of class 'bootmex'.
}
\usage{
coef.bootmex(x, which = "gpd")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class 'bootmex' as returned by \code{bootmex}. }
  \item{which}{ Which set of parameters to return: the parameters from
the generalized Pareto distribution (the default), or the parameters in the 
dependence structure (\code{which="dependence"}). Note that the means
of the bootstrap estimates for the dependence structure parameters are only
meaningful if either all, or no, runs resulted in \code{c} and \code{d}
being estimated. }
}
\details{
The means of the parameter estimates are returned. The procedure estimates
\code{c} and \code{d} in the dependence structure only if \code{a = 0} and
\code{b} is negative. So the parameters have mixtures distributions and
the means are not helpful.
}
\value{
Matrix containing the parameter estimates.
}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth }
\seealso{ \code{\link{migpd}} , \code{\link{mexDependence} },
	\code{\link{bootmex}}, \code{\link{print.mexPrediction}},
	\code{\link{summary.mexPrediction}}, \code{plot.mexPrediction}
}
\examples{
data( winter )
mygpd <- migpd( winter, qu=.7, penalty="none" )
myboot <- bootmex( mygpd , wh="NO", gqu=.7 )
coef( myboot )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
