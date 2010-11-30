\name{migpd}
\alias{migpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit multiple independent generalized Pareto models }
\description{
	Fit multiple independent generalized Pareto models as the first step
  of conditional multivariate extreme values modelling following the
  approach of Heffernan and Tawn, 2004.
}
\usage{
migpd(data, th, qu, penalty = "gaussian", maxit = 10000, verbose = FALSE, gaussianParameters = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A matrix or data.frame, each column of which is to be modelled. }
  \item{th}{ Thresholds above which to fit the models. Only one of
				\code{th} and \code{qu} should be supplied. }
  \item{qu}{ Quantiles above which to fit the models. Only one of
				\code{th} and \code{qu} should be supplied. }
  \item{penalty}{ How the likelihood should be penalized. Defaults to
		'gaussian'. Alternatives are 'none' (resulting in maximum
		likelihood estimation) and 'jeffreys'. }
  \item{maxit}{ The maximum number of iterations to be used by the
		optimizer.}
  \item{verbose}{ Controls whether or not the function print to screen
		every time it fits a model. Defaults to FALSE.}
  \item{gaussianParameters}{ A named list, each element of which contains
			two components: the first should be a vector of length 2
			corresponding to the location of the Gaussian distribution;
			the second should be 2 by 2 matrix corresponding to the 
			covariance of the distribution. This is only used if
			\code{penalty = 'gaussian'}. The names should match the
			names of the columns of \code{data}. If not provided, 
			it defaults to being centred at zero, with variance
			10000 for log(sigma) and 0.25 for xi. See the details section.}
}
\details{
  The parameters in the generalized Pareto distribution are
  estimated for each column of the data in turn, independently
  of all other columns.
  
  Maximum likelihood estimation often fails with generalized
  Pareto distributions because of the likelihood becoming flat
  (see, for example, Hosking et al, 1985).
  Therefore the function allows penalized likelihood estimation,
  which is the same as maximum a posteriori estimation from a 
  Bayesian point of view. 
  
  The Jeffreys prior (the default) is discussed by Castellanos
  and Cabras. An alternative is to use quadratic penalization,
  corresponding to using a Gaussian prior. If no genuine
  prior information is available, the following argument can
  be used. If xi = -1, the generalized Pareto distribution 
  corresponds to the uniform distribution, and if xi is 1 or
  greater, the expectation is infinite. Thefore, xi is likely
  to fall in the region (-1, 1). A Gaussian distribution
  centred at zero and with standard deviation 0.5 will have 
  little mass outside of (-1, 1) and so will often be a
  reasonable prior for xi. For log(sigma) a Gaussian
  distribution, centred
  at zero and with standard deviation 100 will often be vague.
  If a Gaussian penalty is specified but no parameters are
  given, the function will assume such a prior.
  
  Note that internally the function works with log(sigma), not
  sigma. The reasons are that quadratic penalization makes more
  sense for log(sigma) than for sigma (because the distribution
  of log(sigma) will be more nearly symmetric), and because
  it was found to stabilize computations.
  
  The associated \code{coef}, \code{print} and \code{summary}
  functions exponentiate the log(sigma) parameter to return results
  on the expected scale. If you are accessesing the parameters
  directly, however, take care to be sure what scale the results
  are on.
  
  The \code{mex} library contains no functions to help select
  appropriate thresholds or perform diagnostic checks. You should
  use one of the other available extreme values packages to do
  those operations.
}
\value{
  An object of class 'migpd'. There are \code{coef}, \code{print}
  and \code{summary} functions available.
}
\references{
M. E. Castellanos and S. Cabras, A default Bayesian procedure for the generalized Pareto distribution, Journal of Statistical Planning and Inference, 137, 473 -- 483, 2007

J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004

J. R. M. Hosking and J. R. Wallis, Parameter and quantile estimation for the genralized Pareto distribution, Technometrics, 29, 339 -- 349, 1987
}
\author{ Harry Southworth }
\seealso{ \code{\link{mexDependence}}, \code{\link{bootmex}},
	\code{\link{predict.mex}}
		  }
\examples{
data( winter )
mygpd <- migpd( winter, qu=.7, penalty = "none" )
mygpd

}
\keyword{ models }
\keyword{ multivariate }

