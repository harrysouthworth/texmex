\name{mexPrediction}
\alias{mexPrediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prediction for a conditional multivariate exteme values model }
\description{
  Prediction for a conditional multivariate exteme values model following
the approach of Heffernan and Tawn, 2004.
}
\usage{
mexPrediction(migpd, boot, pqu = 0.99, nsim = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{migpd}{ An object of class 'migpd' as returned by function \code{migpd}. }
  \item{boot}{ An object of class 'bootmex' as returned by function \code{bootmex}. }
  \item{pqu}{ The quantile of the conditioning variable over which it will be
thresholded. Defaults to \code{pqu = .99 }. }
  \item{nsim}{ The number of simulated observations to be generated for each
set of bootstrap parameter estimates. }
}
\details{
	For each set of parameter estimates and normalized conditional Gumbel
data, simulated values of the dependent variables are created, given that the
conditioning variable is above its 100\code{pqu} quantile.
}
\value{
An object of class 'mexPrediction'. There are \code{print}, \code{summary} and
\code{plot} functions available. The simulated values based on the actual data
and point estimates appear in component \code{object$data$simulated}. The
simulated data from the bootstrap estimates appears in \code{object$replicates}.
}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth }
\note{ Many thanks to Janet Heffernan for answering lots of questions. }
\seealso{ \code{\link{migpd}} , \code{\link{mexDependence} },
	\code{\link{bootmex}}, \code{\link{print.mexPrediction}},
	\code{\link{summary.mexPrediction}}, \code{plot.mexPrediction}
}
\examples{
data( winter )
mygpd <- migpd( winter, qu=.7 )
myboot <- bootmex( mygpd , which = "NO", gqu=.7 )
mypred <- mexPrediction( mygpd, myboot , pqu=.95 )
summary( mypred )
plot( mypred )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
