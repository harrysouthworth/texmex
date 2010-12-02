\name{plot.mexDependence}
\alias{plot.mexDependence}
\title{ Diagnostic plots for the estimated dependence model in a multivariate conditional dependence extreme value model. }
\description{
Diagnostic plots for dependence model from an object of class \code{mexDependence}.
}
\usage{
plot.mexDependence(x, quantiles=seq(0.1,by=0.2,len=5),col="grey",...)
}
\arguments{
  \item{x}{ An object of class \code{mexDependence}. }
  \item{quantiles}{ A vector of quantiles taking values between 0 and 1 specifying 
  the quantiles of the conditional distributions which will be plotted.}
  \item{col}{ The color for points on scatterplots of residuals and original data. }
}
\details{
    Diagnostic plots for the fitted dependence model described by Heffernan and Tawn, 2004.
    
    Scatterplots of the residuals Z from the fitted model of Heffernan and Tawn (2004) are 
    plotted against the quantile of the conditioning variable, with a lowess curve showing the local 
    mean of these points.  Any trend in the location of these variables with the conditioning variable
    indicates a violation of the model assumption that the residuals Z are indpenendent of the conditioning
    variable.
    
    The absolute value of Z-mean(Z) is also plotted, again with the lowess curve showing 
    the local mean of these points.  These plots are intended to highlight any trend between the 
    variability of the residuals Z and the conditioning variable.
    
    The final plots show the fitted quantiles (specified by \code{quantiles}) of the conditional distribution 
    of each variable given the conditioning variable.  A model that fits well will have good agreement between the
    distribution of the raw data (shown by the scatter plot) and the fitted quantiles.
}
\value{NULL}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach for multivariate
extreme values, Journal of the Royal Statistical Society (B), 66, 497 -- 546, 2004. }
\author{ Harry Southworth }
\seealso{ \code{\link{migpd}}, \code{\link{bootmex}}, \code{\link{predict.mex}} }
\examples{
w <- migpd(winter,mqu=0.7)
wmex <- mexDependence(w)
par(mfrow=c(2,2))
plot(wmex)
}
\keyword{ models }
\keyword{ multivariate }
