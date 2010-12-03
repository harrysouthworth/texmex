\name{mex}
\alias{mex}
Conditional multivariate extreme values modelling
\title{
%%  ~~function to do ... ~~
}
\description{
Fit the conditional multivariate extreme value model of Heffernan and Tawn
}
\usage{
mex(data, which, th, qu, gth, gqu, penalty = "gaussian", maxit = 10000, trace = 0, verbose = FALSE, priorParameters = NULL)
\method{print}{mex}(x, ...)
\method{summary}{mex}(object, ...)
\method{plot}{mex}{x, ...}
\method{coef}{mex}{object, ...}
}

\arguments{
  \item{data}{A numeric matrix or data.frame, the columns of which are to be modelled.}
  \item{which}{The variable on which to condition.}
  \item{mth}{The threshold above which to fit generalized Pareto distributions.
            If this is a vector of length 1, the same threshold will be used for
            each variable. Otherwise, it should be a vector whose length is
            equal to the number of columns in \code{data}.}
  \item{mqu}{As an alternative to specifying the thresholds via th, you can
            specify the quantile above which to fit generalized Pareto
            distributions. If this is a vector of length 1, the same quantile
            will be used for each variable. Otherwise, it should be a vector
            whose length is equal to the number of columns in data.}
  \item{dth}{The threshold to use on the transformed scale; that is, the
            threshold to use when estimating the dependence structure. The same
            threshold will be used for each variable.}
  \item{dqu}{As an alternative to \code{dth}, specify the quantile at which to
            threshold the transformed data when estimating the dependence
            parameters.}
  \item{penalty}{How to penalize the likelihood when estimating the marginal
            generalized Pareto distributions. Defaults to ``gaussian''. See the
            help file for \code{\link{gpd}} for more information.}
  \item{maxit}{The maximum number of iterations to be used by the optimizer.
            defaults to \code{maxit = 1000}.}
  \item{trace}{Whether or not to inform the user of the progress of the
            optimizer. Defaults to 0, indicating no trace.}
  \item{verbose}{Whether or not to keep the user informed of progress. Defaults
            to \code{verbose = FALSE}.}
  \item{priorParameters}{A named list, each element of which contains two
            components: the first should be a vector of length 2 corresponding
            to the location of the Gaussian distribution; the second should be 2
            by 2 matrix corresponding to the covariance of the distribution.
            This is only used if \code{penalty = 'gaussian'}. The names should match
            the names of the columns of \code{data}. If not provided, it defaults to
            being centred at zero, with variance 10000 for log(sigma) and 0.25
            for xi. See the details section.}
}
\details{The parameters in the generalized Pareto distribution are estimated for
         each column of the data in turn, independently of all other columns.
         Then, the conditional multivariate approach of Heffernan and Tawn is
         used to model the dependence between variables.}
\value{
  \item{margins}{An object of class \code{\link{migpd}}.}
  \item{dependence}{An object of class \code{\link{mexDependence}}.}
}
\references{J. E. Heffernan and J. A. Tawn, A conditional approach for
            multivariate extreme values, Journal of the Royal Statistical
            Society B, 66, 497 – 546, 2004}
\author{
Harry Southworth, Janet E. Heffernan
}

\seealso{\code{\link{migpd}}, \code{\link{mexDependence}}, \code{\link{bootmex}}}
\examples{
w <- mex(winter, qu=.7)
w
par(mfcol=c(4, 2))
plot(w)
}
\keyword{ models }
\keyword{ multivariate }
