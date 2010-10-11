\name{plot.cmvxPrediction}
\alias{plot.cmvxPrediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot a 'cmvxPrediction' object. }
\description{
Plot observed and simulated data for a conditional multivariate extreme
values model following Heffernan and Tawn, 2004.
}
\usage{
plot.cmvxPrediction(x, pch = c(1, 3), col = c(2, 8), cex = c(1, 1), ask = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class 'cmvxPrediction' as returned by \code{cmvxPrediction} }
  \item{pch}{ A vector of length 2 giving the plotting characters for the observed
and simulated data respectively. }
  \item{col}{ A vector of length 2 giving the plotting colours for the observed
and simulated data respectively. }
  \item{cex}{ A vector of length 2 giving the plotting character expansion for the observed
and simulated data respectively. }
  \item{ask}{ Should the graphics device ask before clearing the page? By
default, the function checks to see if the device has enough rows and columns
to fit all plots on, and only asks if it does not.}
}
\details{
For each dependent variable, the observed and simulated data are plotted versus
the conditioning variable. A vertical line on the graph indicates the level
at which the conditioning variable was thresholded (in the call to \code{cmvxPrediction}.
Observed and simulated data are differentiated by plotting character and colour and
(possibly) character expansion.
}
\value{
  \code{NULL}
}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth }
\note{ Many thanks to Janet Heffernan for answering lots of questions. }
\seealso{ \code{\link{migpd}} , \code{link{cmvxDependence} },
	\code{\link{cmvxBoot}}, \code{\link{cmvxPrediction}},
	\code{\link{print.cmvxPrediction}},
	\code{\link{summary.cmvxPrediction}}
}
\examples{
data( winter )
mygpd <- migpd( winter, qu=.7 )
myboot <- cmvxBoot( mygpd , which = "NO", gqu=.7 )
mypred <- cmvxPrediction( mygpd, myboot , pqu=.95 )
summary( mypred )
plot( mypred , cex = c( 1, .5 ), col=c( "grey", "blue" ))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
