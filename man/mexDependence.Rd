\name{mexDependence}
\alias{mexDependence}
\alias{print.mexDependence}
\alias{show.mexDependence}
\alias{plot.mexDependence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate the dependence parameters in a conditional multivariate extreme values model }
\description{
  Estimate the dependence parameters in a conditional multivariate extreme
  values model using the approach of Heffernan and Tawn, 2004.
  Many thanks to Janet Heffernan for answering lots of questions.
}
\usage{
mexDependence(x, which, dth, dqu)
\method{print}{mexDependence}(x, ...)
\method{show}{mexDependence}(x, ...)
\method{plot}{mexDependence}(x, quantiles=seq(0.1, by=0.2, len=5), col="grey", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class 'migpd' as returned by \code{migpd}.}
  \item{which}{ The name of the variable on which to condition. This is
		the name of a column of the data that was passed into.
		\code{migpd}}
  \item{dth}{ A single number specifying a threshold on the transformed
			  scale; i.e. the threshold to use when estimating the dependence
              structure. The model will be conditional on the \code{which}
			  variable being over this value. Only one of \code{dth}
			  or \code{dqu} should be given.}
  \item{dqu}{ A single number specifying a quantile on the transformed
			  scale; i.e. the threshold to use when estimating the dependence
              structure. The model will be conditional on the \code{which}
			  variable being over this quantile. Only one of \code{dth}
			  or \code{dqu} should be given.}
  \item{quantiles}{A vector of quantiles taking values between 0 and 1 specifying 
  the quantiles of the conditional distributions which will be plotted.}
  \item{col}{The color for points on scatterplots of residuals and original data.}
  \item{...}{Further arguments to be passed to methods.}
}
\details{
  The parameters a, b, c and d in the dependence structure described
  by Heffernan and Tawn (2004) are estimated. Estimation takes place
  in two steps: first, a and b are estimated; then, if a=0 and b is 
  negative, parameters c and d are estimated. Otherwise c and d will
  be fixed at zero.

  The \code{plot} method produces diagnostic plots for the fitted dependence model described by Heffernan and Tawn, 2004.
    
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
\value{
  An object of class \code{mexDependence}. There are \code{print} and \code{plot}
  functions available.
}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth, Janet E. Heffernan }

\seealso{ \code{\link{migpd}}, \code{\link{bootmex}}, \code{\link{predict.mex}}
		}
\examples{
data(winter)
mygpd <- migpd(winter , mqu=.7, penalty="none")
mexDependence(mygpd , which = "NO", dqu=.7)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }

