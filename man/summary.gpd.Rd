\name{summary.gpd}
\alias{summary.gpd}
\alias{print.gpd}
\alias{print.summary.gpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summarize a \code{gpd} object. }
\description{
 Summaries of an object of class \code{gpd}.
}
\usage{
summary.gpd(x, nsim = 1000, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class \code{gpd}. }
  \item{nsim}{ Number of simulated datasets to use to compute a simulated envelope. }
  \item{alpha}{ The size of the simulation envelope - i.e. a (1 - alpha)\% pointwise
        simulated confidence interval is created. }
}
\details{
  A simulation envelope is produced the data, based on quantiles of the fitted
  model. This represents a pointwise (1 - alpha)\%   simulated confidence interval.
  Since the ordered observations will be correlated, if any observation is outside
  the envelope, it is likely that a chain of observations will be outside the
  envelope. Therefore, if the number outside the envelope is a little more than
  alpha\%, that does not immediately imply a serious shortcoming of the fitted model.
}
\value{
  \item{model}{The \code{gpd} object.}
  \item{envelope }{The simulated envelope.}
  \item{nsim}{The number of simulations.}
  \item{alpha}{The size of the envelope.}
  ...
}

\author{ Harry Southworth }
\seealso{ \code{\link{gpd}}, \code{\link{plot.gpd}} }
\examples{
x <- rnorm(1000)
mod <- gpd(x, qu = 0.7)
mod
summary(mod)
}
\keyword{ models }

