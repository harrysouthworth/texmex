\name{chi}
\alias{chi}
\alias{print.chi}
\alias{summary.chi}
\alias{plot.chi}

\title{ Measures of extremal dependence }
\description{
  Compute measures of extremal dependence for 2 variables.
}
\usage{
chi(data, nq = 100, qlim = NULL, alpha = 0.05, trunc = TRUE)
}

\arguments{
  \item{data}{ A matrix containing 2 numeric columns. }
  \item{nq}{ The number of quantiles at which to evaluate the
             dependence measures. }
  \item{qlim}{ The minimum and maximum quantiles at which to do
               the evaluation. }
  \item{alpha}{ The size of the confidence interval to be
                used. Defaults to \code{alpha = 0.05}. }
  \item{trunc}{ Logical flag indicating whether the estimates
                should be truncated at their theoretical bounds.
                Defaults to \code{trunc = TRUE}. }
}
\details{
  Computes the values chi and chi-bar described by Coles,
  Heffernan and Tawn. The code is a fairly simple reorganization
  of code written by Janet E. Heffernan and Alec Stephenson and
  which appears in
  the \code{chiplot} function in the \code{evd} package.
}
\value{
  An object of class \code{chi} containing the following.

  \item{chi }{Values of chi and their esitmated upper and lower
              confidence limits.}
  \item{chibar }{Values of chibar and their esitmated upper and lower
              confidence limits.}
  \item{quantile}{The quantiles at which chi and chi-bar were
                  evaluated.}
  \item{chiulb, chibarulb}{Uppwer and lower bounds for chi and chi-bar.}
}
\references{ S. Coles, J. E. Heffernan and J. A. Tawn, Dependence
measures for extreme values analyses, Extremes, 2, 339 -- 365, 1999 }
\author{ Janet E. Heffernan, Alec Stephenson, Harry Southworth }

\seealso{ \code{\link{MCS}} }
\examples{

D <- liver[liver$dose == "D",]
chiD <- chi(D[, 5:6])
plot(chiD)
}
\keyword{ multivariate }

