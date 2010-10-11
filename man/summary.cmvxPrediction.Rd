\name{summary.cmvxPrediction}
\alias{summary.cmvxPrediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summarize a 'cmvxPrediction' object. }
\description{
  Summarize a 'cmvxPrediction' object. Expected values and probabilities
of threshold excesses are computed.
}
\usage{
summary.cmvxPrediction(x, th, probs = c(0.05, 0.5, 0.95))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class 'cmvxPrediction', as returned by function
\code{cmvxPrediction}. }
  \item{th}{ Thresholds for the dependent variables. The simulated probability
of excesses over these thresholds will be computed. If omitted, the thresholds
used in the marginal generalized Pareto distributions will be used (i.e. those
that were used by \code{migpd}). }
  \item{probs}{ Quantiles of the distribution of the simulated dependent varialbes.
Defaults to \code{prob = c( .05, .5, .95 )} representinga an empirical bootstrap
median and 90\% confidence interval. }
}

\value{
An object of class 'summary.cmvxPrediction'. There is a \code{print} function
available.}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth }
\seealso{ \code{\link{cmvxPrediction}}, \code{\link{plot.cmvxPrediction}} }
\examples{
data( winter )
mygpd <- migpd( winter, qu = .7, penalty="none" )
myboot <- cmvxBoot( mygpd, gqu=.7, which = "NO" )
mypred <- cmvxPrediction( mygpd, myboot , pqu=.95 )
summary( mypred , probs = c( .025, .5, .975 ) )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
