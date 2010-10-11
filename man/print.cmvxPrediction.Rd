\name{print.cmvxPrediction}
\alias{print.cmvxPrediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print an object of class 'cmvxPrediction' }
\description{
 Print expected values of dependent variables for a conditional 
multivariate extreme values model following the approach of
Heffernan and Tawn, 2004.
}
\usage{
print.cmvxPrediction(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class 'cmvxPrediction'. }
}
\details{
The expected values of the dependent variables, given that the conditioning
variable was above its 100\code{pqu} threshold (specifed in the call to
\code{cmvxPrediction}.
}
\value{
  \code{NULL}
}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth }
\seealso{ \code{\link{migpd}} , \code{\link{cmvxDependence} },
	\code{\link{cmvxBoot}}, \code{\link{cmvxPrediction}},
	\code{\link{summary.cmvxPrediction}}, \code{plot.cmvxPrediction}	    }
\examples{
data( winter )
mygpd <- migpd( winter , qu=.7 )
myboot <- cmvxBoot( mygpd, which = "NO", gqu=.7 )
mypred <- cmvxPrediction( mygpd, myboot, pqu=.95 )
mypred
}
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
