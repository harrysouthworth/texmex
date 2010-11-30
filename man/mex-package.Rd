\name{mex-package}
\alias{mex-package}
\alias{mex}
\docType{package}
\title{
Conditional multivariate extreme values modelling.
}
\description{
Conditional multivariate extreme values modelling following the approach
of Heffernan and Tawn, 2004.

}
\details{
\tabular{ll}{
Package: \tab mex\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-12-06\cr
License: \tab To be decided\cr
}
\code{migpd}: Fit multiple, independent generalized Pareto models to the columns of a data set.

\code{mexDependence}: Estimate the parameters in the dependence structure.

\code{bootmex}: Bootstrap estimation for parameters in generalized Pareto models
and in the dependence structure.

\code{mexPrediction}: Predict values of dependent variables given that one
variable is extreme.

}
\author{
Harry Southworth

Maintainer: Harry Southworth <harry.southworht@astrazeneca.com>

}
\references{
J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004
}
\keyword{models}
\keyword{multivariate}
\keyword{ package }
\examples{
data( winter )
mygpd <- migpd( winter, qu=.7, penalty="none" )
myboot <- bootmex( mygpd , wh="NO", gqu=.7 )
mypred <- mexPrediction( mygpd, myboot , nsim=500 )
plot( mypred )
summary( mypred )
}
