\name{mexRangeFit}
\alias{mexRangeFit}
\title{
Estimate dependence parameters in a conditional multivariate extreme values model over a range of thresholds.
}
\description{
Diagnostic tool to aid the choice of threshold to be used for the estimation of the dependence parameters in the conditional multivariate extreme values model of Heffernan and Tawn, 2004.  
}
\usage{
mexRangeFit(x, which, quantiles = seq(0.5, 0.9, length = 9), R = 10, nPass=3, trace=10, col = "red", bootcol = "grey", margins = "laplace", constrain = TRUE, v = 10, ...)
}
\arguments{
  \item{x}{An object of class \code{\link{mex}} or \code{\link{migpd}}.}
  \item{which}{The variable on which to condition.}
  \item{quantiles}{A numeric vector specifying the quantiles of the marginal distribution of the conditioning variable at which to fit the dependence model.
}
  \item{R}{ The number of bootstrap runs to perform at each threshold. Defaults to \code{R}=10.}
  \item{nPass}{Argument passed to function \code{\link{bootmex}}.}
  \item{trace}{Argument passed to function \code{\link{bootmex}}.}
  \item{col}{
Plotting colour for the point estimates of the dependence parameters.}
  \item{bootcol}{
Plotting colour for the bootstrap estimates of the dependence parameters.}
  \item{margins}{Argument passed to function \code{\link{mexDependence}}.}
  \item{constrain}{Argument passed to function \code{\link{mexDependence}}.}
  \item{v}{Argument passed to function \code{\link{mexDependence}}.}
  \item{\dots}{
Further graphical parameters may be passed, which will be used for plotting.
}
}
\details{
Dependence model parameters are estimated using a range of threshold values.  The sampling variability of these estimates is characterised using the bootstrap.  Point estimates and bootstrap estimates are finally plotted over the range of thresholds.  Choice of threshold should be made such that the point estimates at the chosen threshold and beyond are constant, up to sampling variation.
}
\value{
NULL.
}
\references{
J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004
}
\author{
Harry Southworth, Janet E. Heffernan 
}


\seealso{
\code{\link{mexDependence}}, \code{\link{bootmex}}
}

\examples{
  w <- migpd(winter, mqu=.7)
  w
  par(mfrow=c(4,2))
  mexRangeFit(w,which=1,main="Winter data, Heffernan and Tawn 2004",cex=0.5)
  }

\keyword{ models }
\keyword{ multivariate }

