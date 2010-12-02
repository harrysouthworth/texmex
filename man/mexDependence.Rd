\name{mexDependence}
\alias{mexDependence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate the dependence parameters in a conditional multivariate extreme values model }
\description{
  Estimate the dependence parameters in a conditional multivariate extreme
  values model using the approach of Heffernan and Tawn, 2004.
  Many thanks to Janet Heffernan for answering lots of questions.
}
\usage{
mexDependence(x, which, dth, dqu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class 'migpd' as returned by \code{migpd}. }
  \item{which}{ The name of the variable on which to condition. This is
		the name of a column of the data that was passed into.
		\code{migpd}}
  \item{dth}{ A single number specifying a threshold on the transformed
			  scale; i.e. the threshold to use when estimating the dependence
              structure. The model will be conditional on the \code{which}
			  variable being over this value. Only one of \code{dth}
			  or \code{dqu} should be given. }
  \item{dqu}{ A single number specifying a quantile on the transformed
			  scale; i.e. the threshold to use when estimating the dependence
              structure. The model will be conditional on the \code{which}
			  variable being over this quantile. Only one of \code{dth}
			  or \code{dqu} should be given. }
}
\details{
  The parameters a, b, c and d in the dependence structure described
  by Heffernan and Tawn (2004) are estimated. Estimation takes place
  in two steps: first, a and b are estimated; then, if a=0 and b is 
  negative, parameters c and d are estimated. Otherwise c and d will
  be fixed at zero.
}
\value{
  An object of class \code{mexDependence}. There is a \code{print}
  function available.
}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth }
\note{ Many thanks to Janet Heffernan for answering lots of questions }
\seealso{ \code{\link{migpd}}, \code{\link{bootmex}}, \code{\link{predict.mex}}
		   }
\examples{
data( winter )
mygpd <- migpd(winter , mqu=.7, penalty="none")
mexDependence(mygpd , which = "NO", dqu=.7)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }

