\name{rl}
\alias{rl}
\alias{rl.gpd}
\alias{rl.bgpd}
\alias{rl.bootgpd}
\alias{plot.rl.gpd}
\alias{plot.rl.bgpd}
\alias{plot.rl.bootgpd}
\alias{print.rl.gpd}
\alias{print.rl.bgpd}
\alias{print.rl.bootgpd}
\alias{summary.rl.gpd}
\alias{summary.rl.bgpd}
\alias{summary.rl.bootgpd}

\title{
Return levels
}
\description{
Computation of return levels and confidence intervals for generalized Pareto distributions.
}
\usage{
rl(object, M = 1000, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, 
   alpha = 0.050, unique. = TRUE, ...)
   
\method{rl}{gpd}(object, M = 1000, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, 
   alpha = 0.050, unique. = TRUE, ...)
\method{rl}{bgpd}(object, M = 1000, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, 
   alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL, ...)
\method{rl}{bootgpd}(object, M = 1000, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, 
   alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL, ...)

\method{plot}{rl.gpd}(x, xlab, ylab, main, pch= 1, ptcol =2 , cex=.75, linecol = 4 , 
     cicol = 0, polycol = 15, smooth = TRUE, sameAxes=TRUE, type="median", ... )
\method{plot}{rl.bgpd}(x, xlab, ylab, main, pch= 1, ptcol =2 , cex=.75, linecol = 4 , 
     cicol = 0, polycol = 15, smooth = TRUE, sameAxes=TRUE, type="median", ... )
\method{plot}{rl.bootgpd}(x, xlab, ylab, main, pch= 1, ptcol =2 , cex=.75, linecol = 4 , 
     cicol = 0, polycol = 15, smooth = TRUE, sameAxes=TRUE, type="median", ... )

\method{print}{rl.gpd}(x, digits=3, ...)
\method{print}{rl.bgpd}(x, digits=3, ...)
\method{print}{rl.bootgpd}(x, digits=3, ...)

\method{summary}{rl.gpd}(object, digits=3, ...)
\method{summary}{rl.bgpd}(object, digits=3, ...)
\method{summary}{rl.bootgpd}(object, digits=3, ...)
}
\arguments{

\item{object}{An object of class \code{gpd}, \code{bgpd} or \code{bootgpd}.}
\item{M}{The M-observation return level is computed by the function. Defaults to \code{M = 1000}.}
\item{newdata}{Data from which to calculate the return level. If not provided, the original data
               used to fit the model is used.  Column names must match those of original data matrix used for model fitting.}
\item{se.fit}{Whether or not to return the standard error of the predicted value.
              Defaults to \code{se.fit = FALSE}.}
\item{ci.fit}{Whether or not to return a confidence interval for the predicted
value. Defaults to \code{ci.fit = FALSE}. For objects of class \code{gpd}, if set to \code{TRUE}
then the confidence interval is a simple symmetric confidence interval
based on the estimated approximate standard error. For the \code{bgpd} and \code{bootgpd} methods, the confidence
interval represents quantiles of the simulated distribution of the
parameters.
}
  \item{alpha}{
If \code{ci.fit = TRUE}, a (1 - alpha)\% confidence interval is returned.
Defaults to \code{alpha = 0.050}.
}
  \item{unique.}{
If \code{unique. = TRUE}, predictions for only the unique values of
the linear predictors are returned, rather than for every row of the original dataframe or of
\code{newdata} if this latter is specified. Defaults to \code{unique. = TRUE}.
}
\item{all}{
For the \code{bgpd} and \code{bootgpd} methods, if \code{all = TRUE}, the
predictions are returned for every simulated parameter vector. Otherwise,
only a summary of the posterior/bootstrap distribution is returned.
Defaults to \code{all = FALSE}.
}
\item{sumfun}{
For the \code{bgpd} and \code{bootgpd} methods, a summary function
can be passed in. If \code{sumfun = FALSE}, the default, the
summary function used returns the estimated mean and median, and quantiles
implied by \code{alpha}.
    }
\item{type}{ For calls to plot methods for objects of class \code{rl.bgpd} or \code{rl.bootgpd}, specifies whether to use the sample mean (\code{type="mean"}) or median (\code{type="median"}) estimate of the return levels.}
\item{x}{Object passed to plot and print methods.}
\item{xlab, ylab, main, pch, ptcol, cex, linecol, cicol, polycol, smooth, sameAxes}{Further arguments to plot methods.}
\item{digits}{Number of digits to show when printing output.}
\item{...}{Further arguments to be passed to methods.}
}

\details{
The M-observation return level is defined as the value that is expected to be exceeded only
once every M observations. Thus, it is an estimate of a high quantile of the fitted distribution.
In models fit by the \code{gpd} family of functions, only a fraction of the data is actually
included in the model; the fitted GPD is a conditional model, conditioning on the threshold having been exceeded. This consideration is taken into account by \code{rl} which calculates unconditional return levels from the entire distribution of observations above and below the GPD fitting threshold.
}
\examples{
mod <- gpd(rain, qu=.8) # daily rainfall observations
rl(mod, M=100*365) # 100-year return level
}
\keyword{datasets}
