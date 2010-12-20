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
\method{print}{chi}(x, ...)
\method{summary}{chi}(object, digits=3, ...)
\method{plot}{chi}(x, which=1:2, lty=1, cilty=2, col=1, spcases=FALSE,
                   cicol=1, xlim=c(0, 1), ylim1=c(-1, 1), ylim2=c(-1, 1),
                   main1="Chi", main2="Chi Bar", xlab="Quantile",
                   ylab1=expression(chi),
                   ylab2=expression(bar(chi)),
                   ask=nb.fig < length(which) && dev.interactive(), ...)
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
  \item{x, object}{An object of class \code{chi}.}
  \item{digits}{ Number of digits for printing. }
  \item{which}{Which plot to do: 1 results in a plot of chi, 2 in a 
               plot of chi bar. Defaults to \code{which = 1:2}.}
  \item{lty, cilty, col, cicol}{Line types and colours for the the
               estimated quantities and thier confidence intervals.}
  \item{xlim, ylim1, ylim2}{Limits for the axes.}
  \item{main1, main2}{Main titles for the plots.}
  \item{xlab, ylab1, ylab2}{Axis labels for the plots.}
  \item{spcases}{Whether or not to plot special cases. Defaults to \code{FALSE}.}
  \item{ask}{Whether or not to ask before reusing the graphics device.}
  \item{...}{Further arguments to be passed to methods.}
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

\note{When the data contain ties, the values of chi and chibar are calculated by 
assigning distinct ranks to tied values using the \code{rank} function with 
argument \code{ties.method = "first"}.  This results in the values of chi 
and chibar being sensitive to the order in which the tied values appear 
in the data.}

\seealso{ \code{\link{MCS}}, \code{\link{rank}} }
\examples{

D <- liver[liver$dose == "D",]
chiD <- chi(D[, 5:6])
plot(chiD)
}
\keyword{ multivariate }

