\name{bootmex}
\alias{bootmex}
\alias{print.bootmex}
\alias{plot.bootmex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bootstrap a conditional multivariate extreme values model }
\description{
  Bootstrap a conditional multivariate extreme values model following
the method of Heffernan and Tawn, 2004. I'm grateful to Jan Heffernan
for advice. The methods for the object returned by \code{bootmex}
are currently flawed!
}
\usage{
bootmex(x, which, B = 100, gth, gqu, nPass=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class 'migpd' as returned by function \code{migpd}. }
  \item{which}{ The name of the variable on which the conditioning is done. This should
be the name of one of the columns in the data set that was passed to \code{migpd}. }
  \item{B}{ The number of bootstrap runs to perform. Defaults to B=100. }
  \item{gth}{ The threshold on the Gumbel scale to apply to the specified conditioning
variable. Only one of \code{gth} or \code{gqu} should be provided. }
  \item{gqu}{ The quantile on the Gumbel scale to threshold to the specified conditioning
variable at. Only one of \code{gth} or \code{gqu} should be provided. }
  \item{nPass}{An integer. Sometimes, particularly with small samples,
          the estimation process fails with
          some bootstrap samples. The function has a look to see which runs
          fail and takes additional bootstrap samples in an attempt to get
          parameter estimates. By default, it has nPass=3 attempts at this
          before giving up.}
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
}
\value{
 An object of class 'bootmex'. A print function is available.
}

\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth }
\note{ Many thanks to Janet Heffernan for answering lots of questions. }
\seealso{ \code{\link{migpd}} , \code{link{mexDependence} },
	\code{\link{bootmex}}, \code{\link{mexPrediction}},
		    }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data( winter )
mygpd <- migpd( winter , qu = .7 )
myboot <- bootmex( mygpd, which = "NO", gqu=.7 )
myboot
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
