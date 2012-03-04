\name{chi}
\alias{chi}
\alias{summary.chi}
\alias{plot.chi}

\title{ Measures of extremal dependence }
\description{
  Compute measures of extremal dependence for 2 variables.
}
\usage{
chi(data, nq = 100, qlim = NULL, alpha = 0.05, trunc = TRUE)

\method{summary}{chi}(object, digits=3, ...)

\method{plot}{chi}(x, show=c("Chi"=TRUE,"ChiBar"=TRUE), lty=1, cilty=2, col=1, spcases=TRUE,
    cicol=1, xlim=c(0, 1), ylimChi = c(-1, 1), ylimChiBar = c(-1, 1),
    mainChi = "Chi", mainChiBar = "Chi Bar", xlab = "Quantile", 
    ylabChi    = expression(chi(u)),#"Chi(u)"
    ylabChiBar = expression(bar(chi)(u)), #"Chi Bar(u)"
    ask, ...)
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
  \item{show}{Logical, of length 2, names "Chi" and "ChiBar". 
              Defaults to \code{c("Chi" = TRUE, "ChiBar" = TRUE)}.}
  \item{lty, cilty, col, cicol}{Line types and colours for the the
               estimated quantities and their confidence intervals.}
  \item{xlim, ylimChi, ylimChiBar}{Limits for the axes.}
  \item{mainChi, mainChiBar}{Main titles for the plots.}
  \item{xlab, ylabChi, ylabChiBar}{Axis labels for the plots.}
  \item{spcases}{Whether or not to plot special cases of perfect (positive and negative) dependence and indpenendence. Defaults to \code{FALSE}.}
  \item{ask}{Whether or not to ask before reusing the graphics device.}
  \item{...}{Further arguments to be passed to methods.}
}
\details{
  Computes the functions chi and chi-bar described by Coles,
  Heffernan and Tawn. The limiting values of these functions as the 
  \code{Quantile} approaches 1 give an empirical measure of the type 
  and strength of tail dependendce exhibited by the data. 
  
  A limiting value of ChiBar equal to 1 indicates Asymptotic Dependence, 
  in which case the limiting value of Chi gives a measure of the strength 
  of dependence in this class.  A limiting value of ChiBar of less than 1 
  indicates Asymptotic Independence in which case the limiting value of 
  ChiBar gives a measure of the strength of dependence. 
  
  The plot method shows the ChiBar and Chi functions.  In the case of the 
  confidence interval for ChiBar excluding the value 1 for all of the largest 
  values of \code{Quantile}, the plot of the Chi function is shown in grey.
}
\value{
  An object of class \code{chi} containing the following.

  \item{chi }{Values of chi and their esitmated upper and lower
              confidence limits.}
  \item{chibar }{Values of chibar and their esitmated upper and lower
              confidence limits.}
  \item{quantile}{The quantiles at which chi and chi-bar were
                  evaluated.}
  \item{chiulb, chibarulb}{Upper and lower bounds for chi and chi-bar.}
}
\references{ S. Coles, J. E. Heffernan and J. A. Tawn, Dependence
measures for extreme values analyses, Extremes, 2, 339 -- 365, 1999. }
\author{ Janet E. Heffernan, Alec Stephenson, Harry Southworth }

\note{When the data contain ties, the values of chi and chibar are calculated by 
assigning distinct ranks to tied values using the \code{rank} function with 
argument \code{ties.method = "first"}.  This results in the values of chi 
and chibar being sensitive to the order in which the tied values appear 
in the data.

The code is a fairly simple reorganization
  of code written by Janet E. Heffernan and Alec Stephenson and
  which appears in
  the \code{chiplot} function in the \code{evd} package.
  }

\seealso{ \code{\link{MCS}}, \code{\link{rank}} }
\examples{

D <- liver[liver$dose == "D",]
chiD <- chi(D[, 5:6])
par(mfrow=c(1,2))
plot(chiD)

A <- liver[liver$dose == "A",]
chiA <- chi(A[, 5:6])
# here the limiting value of chi bar(u) lies away from one so the chi plot is not relevant and is plotted in grey
plot(chiA) 


}
\keyword{ multivariate }

