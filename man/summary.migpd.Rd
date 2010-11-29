\name{summary.migpd}
\alias{summary.migpd}
\alias{print.migpd}
\alias{print.summary.migpd}
\alias{coef.migpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summaries of multiple independent fitted Generalized Pareto Distributions }
\description{
Summaries of an object of class \code{migpd}.
}
\usage{
summary(x)
print(x)
coef(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class \code{migpd}. }
}
\details{
    No plot function is provided. The user should use the \code{gpd} function
    and look at its output to establish if the marginal models are good fits.
}
\value{
A matrix containing the marginal point estimates and the thresholds.
}

\author{ Harry Southworth }
\seealso{
    \code{\link{gpd}}, \code{\link{mexDependence}}, \code{\link{migpdCoefs}}
}
\examples{
mod <- migpd(winter, qu=.7)
mod
summary(mod)
coef(mod)
}
\keyword{ models }
\keyword{ multivariate }
