\name{bootmex}
\alias{bootmex}
\alias{print.bootmex}
\alias{show.bootmex}
\alias{plot.bootmex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bootstrap a conditional multivariate extreme values model }
\description{
  Bootstrap a conditional multivariate extreme values model following
the method of Heffernan and Tawn, 2004. I'm grateful to Jan Heffernan
for advice.
}
\usage{
bootmex(x, which, R = 100, dth, dqu, nPass=3, trace=10)
\method{plot}{bootmex}(x, plots = "gpd", main = "", ...)
\method{print}{bootmex}(x, ...)
\method{show}{bootmex}(x, ...)
}
\arguments{
  \item{x}{ An object of class 'migpd' as returned by function \code{migpd}. }
  \item{which}{ The name of the variable on which the conditioning is done. This should
be the name of one of the columns in the data set that was passed to \code{migpd}. }
  \item{R}{ The number of bootstrap runs to perform. Defaults to \code{R}=100. }
  \item{dth}{ The threshold on the Gumbel scale to apply to the specified conditioning
variable. Only one of \code{dth} or \code{dqu} should be provided. If \code{x} has
        class \code{mex}, the value is read from the object and need not be
        provided.}
  \item{dqu}{ The quantile to threshold to the specified conditioning
variable at. Only one of \code{dth} or \code{dqu} should be provided. If \code{x} has
        class \code{mex}, the value is read from the object and need not be
        provided.}
  \item{nPass}{An integer. Sometimes, particularly with small samples,
          the estimation process fails with
          some bootstrap samples. The function has a look to see which runs
          fail and takes additional bootstrap samples in an attempt to get
          parameter estimates. By default, it has nPass=3 attempts at this
          before giving up.}
  \item{trace}{How often the inform the user of progress. Defaults to \code{trace=10}.}
  \item{plots}{ Whether to do diagnostic plots of gpd parameter estimates or
      dependence parameter estimates.  Defaults to "gpd". }
  \item{main}{Title for plots.}
  \item{...}{Further arguments to be passed to methods.}
}
\details{
 Details of the bootstrap method are given by Heffernan and Tawn (2004). Firstly, Gumbel
data are simulated independently. A bootstrap sample of the Gumbel-transformed original
data is then taken and the simulated Gumbel variables are sorted to match the sort order
of this bootstrap sample. This procedure ensures the simulated data has the same dependence
structure as the bootstrap sample.

The simulated, ordered Gumbel data are then transformed to the scale of the original data.
Values beneath the original threshold are passed through the empirical distribution function
and values above that threshold are passed through the appropriate distribution function
of the estimated generalized Pareto distribution. The parameters in the generalized Pareto
distributions are then estimated, the data are transformed to the Gumbel scale using the
orginal threshold, the empirical distribution function and these estimated parameters. The
variables in the dependence structure of these variables are then estimated.

Note that maximum likelihood estimation will often fail for small samples when the generalized
Pareto distribution is being fit. Therefore it will often be useful to use penalized likelihood
estimation. \code{bootmex} does whatever was done in the call to \code{migpd}.

Also note that sometimes (again, usually with small data sets) all of the simulated Gumbel
random numbers will be beneath the threshold for the conditioning variable. Such samples are
abandoned by \code{bootmex} and a new sample is generated. This probably introduces some
bias into the resulting bootstrap distributions.

The \code{plot} method produces histograms of bootstrap gpd parameters (the default)
or scatterplots of dependence parameters with the point estimates for the
original data shown.

By design, there is no \code{coef} method. The bootstrapping is done to account
for uncertainty. It is not obvious that adjusting the parameters for the
mean bias is the correct thing to do.
}
\value{
 An object of class 'bootmex'. Print and plot functions are available.
}

\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth }
\note{ Many thanks to Janet Heffernan for answering lots of questions. }
\seealso{ \code{\link{migpd}} , \code{link{mexDependence} },
	\code{\link{bootmex}}, \code{\link{predict.mex}}.
		    }
\examples{
# Uncomment the following lines to run example - commented out to keep CRAN robots happy
#mygpd <- migpd(winter , mqu = .7)
#myboot <- bootmex(mygpd, which = "NO", dqu=.7)
#myboot
#plot(myboot,plots="gpd")
#plot(myboot,plots="dependence")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
