\name{predict.evmOpt}
\alias{predict.evmOpt}
\alias{predict.evmSim}
\alias{predict.evmBoot}
\alias{linearPredictors}
\alias{linearPredictors.evmOpt}
\alias{linearPredictors.evmSim}
\alias{linearPredictors.evmBoot}
\alias{plot.lp.evmOpt}
\alias{plot.lp.evmSim}
\alias{plot.lp.evmBoot}
\alias{print.lp.evmOpt}
\alias{print.lp.evmSim}
\alias{print.lp.evmBoot}
\alias{summary.lp.evmOpt}
\alias{summary.lp.evmSim}
\alias{summary.lp.evmBoot}
\alias{predict.evmSim}

\title{Predict return levels from extreme value models,
       or obtain the linear predictors.}
\description{Predict return levels from extreme value models, 
       or obtain the linear predictors.}

\usage{
\method{predict}{evmOpt}(object, M = 1000, newdata = NULL,
type = "return level", se.fit = FALSE,
ci.fit = FALSE, alpha = 0.05, unique. = TRUE,...)

\method{predict}{evmSim}(object, M = 1000, newdata = NULL,
type = "return level", se.fit = FALSE,
ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE,
sumfun = NULL,...)

\method{predict}{evmBoot}(object, M = 1000, newdata = NULL,
type = "return level", se.fit = FALSE,
ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE,
sumfun = NULL,...)

linearPredictors(object, newdata = NULL, se.fit = FALSE, ci.fit = FALSE,
                 alpha = 0.050, unique. = TRUE, ...)

\method{linearPredictors}{evmOpt}(object, newdata = NULL, se.fit = FALSE, ci.fit = FALSE,
alpha = 0.05, unique. = TRUE, full.cov = FALSE,...)
\method{linearPredictors}{evmSim}(object, newdata = NULL, se.fit = FALSE, ci.fit = FALSE,
alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL,...)
\method{linearPredictors}{evmBoot}(object, newdata = NULL, se.fit = FALSE,
                                   ci.fit = FALSE, alpha = 0.050, unique. = TRUE,
                                   all = FALSE, sumfun = NULL,...)

\method{print}{lp.evmOpt}(x, digits=3,...)
\method{print}{lp.evmSim}(x, digits=3,...)
\method{print}{lp.evmBoot}(x, digits=3,...)

\method{summary}{lp.evmOpt}(object, digits=3,...)
\method{summary}{lp.evmSim}(object, digits=3,...)
\method{summary}{lp.evmBoot}(object, digits=3,...)

\method{plot}{lp.evmOpt}(x, main=NULL, pch=1, ptcol=2, cex=.75,
                         linecol=4, cicol=1, polycol=15,...)
\method{plot}{lp.evmSim}(x, type="median", ...)
\method{plot}{lp.evmBoot}(x, type="median", ...)
}
\arguments{
  \item{object}{
An object of class \code{evmOpt}, \code{evmSim} or \code{evmBoot}.
}
  \item{newdata}{
The new data that you want to make the prediction for. Defaults in
\code{newdata = NULL} in which case the data used in fitting the model
will be used. Column names must match those of the original data
matrix used for model fitting.}

  \item{type}{For the predict methods, the type of prediction, either
      "return level" (or "rl") or "link" (or "lp"). Defaults to
      \code{type = "return level"}. When a return level is wanted, the
      user can specify the associated return period via the \code{M} 
      argument. If \code{type = "link"} the linear predictor(s) for
      \code{phi} and \code{xi} (or whatever other parameters are
      in your \code{texmexFamily} are returned. 

      For the plot methods for simulation based estimation of 
      underlying distributions i.e. objects derived from "evmSim" and 
      "evmBoot" classes, whether to use the sample median
      \code{type="median"} or mean \code{type="mean"} estimate of the 
      parameter.
}
  \item{se.fit}{
Whether or not to return the standard error of the predicted value.
Defaults to \code{se.fit = FALSE} and is not implemented for
\code{predict.evmSim} or \code{predict.evmBoot}.
}
  \item{ci.fit}{
Whether or not to return a confidence interval for the predicted
value. Defaults to \code{ci.fit = FALSE}. For objects of class
\code{evmOpt}, if set to \code{TRUE}
then the confidence interval is a simple symmetric confidence interval
based on the estimated approximate standard error. For the \code{evmSim}
and \code{evmBoot} methods, the confidence
interval represents quantiles of the simulated distribution of the
parameters.
}
  \item{M}{
The return level: units are number of observations. Defaults to
\code{M = 1000}. If a vector is passed, a list is returned, with items
corresponding to the different values of the vector \code{M}.
}
  \item{alpha}{
If \code{ci.fit = TRUE}, a 100(1 - alpha)\% confidence interval is returned.
Defaults to \code{alpha = 0.050}.
}
  \item{unique.}{
If \code{unique. = TRUE}, predictions for only the unique values of
the linear predictors are returned, rather than for every row of
\code{newdata}. Defaults to \code{unique. = TRUE}.
}
   \item{all}{
For the \code{evmSim} and \code{evmBoot} methods, if \code{all = TRUE}, the
predictions are returned for every simulated parameter vector. Otherwise,
only a summary of the posterior/bootstrap distribution is returned.
Defaults to \code{all = FALSE}.
}
    \item{full.cov}{
Should the full covariance matrix be returned as part of a \code{list}
object. This is used internally and not intended for direct use.
Defaults to \code{full.cov = FALSE}
    }
    \item{sumfun}{
For the \code{evmSim} and \code{evmBoot} methods, a summary function
can be passed in. If \code{sumfun = FALSE}, the default, the
summary function used returns the estimated mean and median, and quantiles
implied by \code{alpha}.
    }
    \item{x}{An object of class \code{lp.evmOpt}, \code{lp.evmSim}
       or \code{lp.evmBoot}, to be passed to methods for these classes.}
    \item{main, pch, ptcol, cex, linecol, cicol, polycol}{Further arguments to plot methods.}
    \item{digits}{Number of digits to show when printing objects.}
    \item{...}{Further arguments to methods.}
}
\details{
By default, return levels predicted from the unique values of the
linear predictors are returned. For the \code{evmBoot} method,
estimates of confidence intervals are simply quantiles of the bootstrap
sample. The \code{evmBoot}
method is just a wrapper for the \code{evmSim} method.
}
\value{
A list with one entry for each value of \code{M}.}
\author{Harry Southworth and Janet E. Heffernan}
\note{
At present, the confidence intervals returned for an object of class
\code{evmOpt} are simple confidence intervals based on assumptions
of normality that are likely to be far from the truth in many cases.
A better approach would be to use profile likelihood, and we intend
to implement this method at a future date.  Alternatively, the
credible intervals returned by using Bayesian estimation and the
predict method for class "evmSim" will tend to give a better
representation of the asymmetry of the estimated intervals around
the parameter point estimates.
}
\keyword{ methods }
