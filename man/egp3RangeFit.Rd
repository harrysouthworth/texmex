\name{egpRangeFit}
\alias{egp3RangeFit}
\alias{print.egp3RangeFit}
\alias{plot.egp3RangeFit}

\title{Estimate the EGP3 distribution power parameter over a range of thresholds}
\description{
  Estimate extended generalized Pareto distribution power parameter over a range of values, using
  maximum (penalized) likelihood.
}
\usage{
egp3RangeFit(data, umin=quantile(data, .05), umax=quantile(data, .95), nint = 10, 
            penalty = "gaussian", priorParameters = NULL, alpha=0.05)
\method{print}{egp3RangeFit}(x, ...)
\method{summary}{egp3RangeFit}(object, ...)
\method{plot}{egp3RangeFit}(x, xlab = "Threshold", ylab = "kappa", main = NULL, addNexcesses=TRUE, ...)
}

\arguments{
  \item{data}{ The data vector to be modelled. }
  \item{umin}{ The minimum threshold above which to estimate the parameters. }
  \item{umax}{ The maximum threshold above which to estimate the parameters. }
  \item{nint}{ The number of thresholds at which to perform the estimation. }
  \item{penalty}{The type of penalty to be used in the maximum penalized
               likelihood estimation. Should be either "gaussian" or 
               "none". Defaults to "gaussian".}
  \item{priorParameters}{Parameters to be used for the penalty function.
                 See the help for \code{\link{evm}} for more informaiton.}
  \item{alpha}{100(1 - alpha)\% confidence intervals will be plotted with the 
               point estimates. Defaults to \code{alpha = 0.05}.}
  \item{x, object}{Arguments to \code{print} and \code{summary} functions.}
  \item{xlab}{ Label for the x-axis. }
  \item{ylab}{ Label for the y-axis. }
  \item{main}{ The main title. }
  \item{addNexcesses}{ Annotate top axis with numbers of threshold excesses
                       arising with the corresponding values of threshold on the
                       bottom axis.}
  \item{\dots}{ Arguments to \code{plot} }
}
\details{
  Papastathopoulos and Tawn present 3 extended versions of the generalized Pareto
  distribution. Using the \code{egp3} texmex family object, the power parameter
  in the EGP3 distribution is estimated on the log scale, a confidence interval is
  calculated and the result is transformed back to the scale of the power parameter
  and returned to the user.
  
  When the power paramer, kappa, is equal to 1, the EPG3 distribution is identical
  to the generalized Pareto distribution. Therefore, the plot of the estimated
  parameter over a range of thresholds provides a diagnostic for threshold selection:
  the lowest value of kappa whose confidence interval includes 1 is suggested as
  the threshold for generalized Pareto modelling.
  
  If lower thresholds are used and the EGP3 distribution itself is used for modelling,
  some care should be taken to ensure the model provides a reasonable degree of
  fit to the data. Limited experience suggests that such models seldom fit well
  and the main value of the EGP3 distribution is as a diagnostic for threshold
  selection as described here.

  Note this function does not extend to assessing model fit when there are
  covariates included in the model.
}

\author{Harry Southworth}
\seealso{ \code{\link{evm}}, \code{\link{gpdRangeFit}}, \code{\link{mrl}} }
\examples{
plot(egp3RangeFit(rain))
}
\keyword{ models }

