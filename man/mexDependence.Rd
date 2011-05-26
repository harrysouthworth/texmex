\name{mexDependence}
\alias{mexDependence}
\alias{print.mexDependence}
\alias{show.mexDependence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate the dependence parameters in a conditional multivariate extreme values model }
\description{
  Estimate the dependence parameters in a conditional multivariate extreme
  values model using the approach of Heffernan and Tawn, 2004.
  Many thanks to Janet Heffernan for answering lots of questions.
}
\usage{
mexDependence(x, which, dth, dqu,margins="laplace", constrain=TRUE, v = 10, maxit=10000, start=c(.01, .01), nGridPlotLik=NULL, marTransform="mixture")
\method{print}{mexDependence}(x, ...)
\method{show}{mexDependence}(x, ...)
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
  \item{margins}{The form of margins to which the data are transformed for carrying out dependence estimation.  Defaults to "laplace", with the alternative option being "gumbel".  The choice of margins has an impact on the interpretation of the fitted dependence parameters.  Under Gumbel margins, the estimated parameters a and b describe only positive dependence, while c and d describe negative dependence in this case.  For Laplace margins, only parameters a and b are estimated as these capture both positive and negative dependence.}
  \item{constrain}{Logical value.  Defaults to \code{constrain=TRUE} although this will subjequently be changed to FALSE if \code{margins="gumbel"} for which constrained estimation is not implemented.  If \code{margins="laplace"} and \code{constrain=TRUE} then the dependence parameter space is constrained to allow only combinations of parameters which give the correct stochastic ordering between (positively and negatively) asymptotically dependent variables and variables which are asymptotically independent.}
  \item{v}{Scalar. Tuning parameter used to carry out constrained estimation of dependence structure under \code{constrain=TRUE}. Takes positive values greater than 1; values between 2 and 10 are recommended.}
  \item{maxit}{The maximum number of iterations to be used by the optimizer.
            defaults to \code{maxit = 1000}.}
  \item{start}{Optional starting value for dependence estimation.  This can be: a vector of length two, with values corresponding to dependence parameters a and b respectively, and in which case \code{start} is used as a starting value for numerical estimation of each of the dependence models to be estimated; a matrix with two rows corresponding to dependence parameters a and b respectively and number of columns equal to the number of dependence models to be estimated (the ordering of the columns will be as in the original data matrix); or a previously estimated object of class \code{'mex'} whose dependence parameter estimates are used as a starting point for estimation.  Note that under \code{constrain=TRUE}, if supplied, \code{start} must lie within the permitted area of the parameter space.} 
  \item{nGridPlotLik}{Optional integer.  Dimension of square grid on which to plot the profile likelihood surface for parameters a and b.  Defaults to \code{nGridPlotLik=NULL} in which case no plotting is carried out.}
  \item{marTransform}{Optional form of transformation to be used for probability integral transform of data from original to Gumbel or Laplace margins.  Takes values \code{marTransform="mixture"} (the default) or \code{marTransform="empirical"}. When \code{marTransform="mixture"}, the rank transform is used below the corresponding GPD fitting threshold used in \code{x}, and the fitted gpd tail model is used above this threshold.  When \code{marTransform="empirical"} only the rank transform is used.}
  \item{...}{Further arguments to be passed to methods.}
}
\details{
  Estimates the extremal dependence structure of the data in \code{x}.  The precise nature of the estimation depends on the value of \code{margins}.  If \code{margins="laplace"} (the default) then dependence parameters a and b are estimated after transformation of the data to Laplace marginal distributions.  These parameters can describe both positive and negative dependence.   If \code{margins="gumbel"} then
  the parameters a, b, c and d in the dependence structure described
  by Heffernan and Tawn (2004) are estimated in the following two steps: first, a and b are estimated; then, if a=0 and b is 
  negative, parameters c and d are estimated (this is the case of negative dependence). Otherwise c and d will
  be fixed at zero (this is the case of positive dependence). 

}
\value{
  An object of class \code{mexDependence}. There is a \code{print} method available.
}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth, Janet E. Heffernan }

\seealso{ \code{\link{migpd}}, \code{\link{bootmex}}, \code{\link{predict.mex}}, \code{\link{plot.mex}}
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

