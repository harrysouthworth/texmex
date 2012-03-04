\name{predict.gpd}
\alias{predict.gpd}
\alias{predict.bgpd}
\alias{predict.bootgpd}
\alias{linearPredictors}
\alias{linearPredictors.gpd}
\alias{linearPredictors.bgpd}
\alias{linearPredictors.bootgpd}
\alias{plot.lp.gpd}
\alias{plot.lp.bgpd}
\alias{plot.lp.bootgpd}
\alias{print.lp.gpd}
\alias{print.lp.bgpd}
\alias{print.lp.bootgpd}
\alias{summary.lp.gpd}
\alias{summary.lp.bgpd}
\alias{summary.lp.bootgpd}
\alias{predict.bgpd}

\title{Predict return levels from Generalized Pareto Distribution models, 
       or obtain the linear predictors.}
\description{Predict return levels from Generalized Pareto Distribution models, 
       or obtain the linear predictors.}

\usage{
\method{predict}{gpd}(object, M = 1000, newdata = NULL, type = "return level", se.fit = FALSE, 
ci.fit = FALSE, alpha = 0.05, unique. = TRUE,...)

\method{predict}{bgpd}(object, M = 1000, newdata = NULL, type = "return level", se.fit = FALSE, 
ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL,...)

\method{predict}{bootgpd}(object, M = 1000, newdata = NULL, type = "return level", se.fit = FALSE, 
ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL,...)


linearPredictors(object, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, 
unique. = TRUE, ...)

\method{linearPredictors}{gpd}(object, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, 
alpha = 0.05, unique. = TRUE, full.cov = FALSE,...)
\method{linearPredictors}{bgpd}(object, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, 
alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL,...)
\method{linearPredictors}{bootgpd}(object, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL,...)

\method{print}{lp.gpd}(x, digits=3,...)
\method{print}{lp.bgpd}(x, digits=3,...)
\method{print}{lp.bootgpd}(x, digits=3,...)

\method{summary}{lp.gpd}(object, digits=3,...)
\method{summary}{lp.bgpd}(object, digits=3,...)
\method{summary}{lp.bootgpd}(object, digits=3,...)

\method{plot}{lp.gpd}(x, main=NULL, pch=1, ptcol=2, cex=.75, linecol=4, cicol=1, polycol=15,...)
\method{plot}{lp.bgpd}(x, type="median", ...)
\method{plot}{lp.bootgpd}(x, type="median", ...)
}
\arguments{
  \item{object}{
An object of class \code{gpd}, \code{bgpd} or \code{bootgpd}.
}
  \item{newdata}{
The new data that you want to make the prediction for. Defaults in
\code{newdata = NULL} in which case the data used in fitting the model
will be used. Column names must match those of original data matrix used for model fitting.}

  \item{type}{
For the predict methods, the type of prediction, either ``return level'' (or ``rl'') or
``link'' (or ``lp''). Defaults to \code{type = "return level"}. When
a return level is wanted, the user can specify the associated return perdiod via the \code{M} argument. If \code{type = ``link''} the linear
predictor(s) for \code{phi} and \code{xi} are returned. 

For the plot methods for simulation based estimation of underlying distributions i.e. objects derived from bgpd and bootgpd classes, whether to use the sample median \code{type="median"} or mean \code{type="mean"} estimate of the parameter.
}
  \item{se.fit}{
Whether or not to return the standard error of the predicted value.
Defaults to \code{se.fit = FALSE} and is not implemented for
\code{predict.bgpd} or \code{predict.bootgpd}.
}
  \item{ci.fit}{
Whether or not to return a confidence interval for the predicted
value. Defaults to \code{ci.fit = FALSE}. For objects of class \code{gpd}, if set to \code{TRUE}
then the confidence interval is a simple symmetric confidence interval
based on the estimated approximate standard error. For the \code{bgpd} and \code{bootgpd} methods, the confidence
interval represents quantiles of the simulated distribution of the
parameters.
}
  \item{M}{
The return level: units are number of observations. Defaults to \code{M = 1000}. If a vector is passed,
a list is returned, with items corresponding to the different values of the vector \code{M}.
}
  \item{alpha}{
If \code{ci.fit = TRUE}, a (1 - alpha)\% confidence interval is returned.
Defaults to \code{alpha = 0.050}.
}
  \item{unique.}{
If \code{unique. = TRUE}, predictions for only the unique values of
the linear predictors are returned, rather than for every row of
\code{newdata}. Defaults to \code{unique. = TRUE}.
}
   \item{all}{
For the \code{bgpd} and \code{bootgpd} methods, if \code{all = TRUE}, the
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
For the \code{bgpd} and \code{bootgpd} methods, a summary function
can be passed in. If \code{sumfun = FALSE}, the default, the
summary function used returns the estimated mean and median, and quantiles
implied by \code{alpha}.
    }
    \item{x}{An object of class \code{lp.gpd}, \code{lp.bgpd} or \code{lp.bootgpd}, to be passed to methods for these classes.}
    \item{main, pch, ptcol, cex, linecol, cicol, polycol}{Further arguments to plot methods.}
    \item{digits}{Number of digits to show when printing objects.}
    \item{...}{Further arguments to methods.}
}
\details{
By default, return levels predicted from the unique values of the
linear predictors are returned. For the \code{bootgpd} method,
estimates of confidence intervals are simply quantiles of the bootstrap sample. The \code{bootgpd}
method is just a wrapper for the \code{bgpd} method.
}
\value{
A list with one entry for each value of \code{M}.}
\author{
Harry Southworth and Janet E. Heffernan
}
\note{
At present, the confidence intervals returned for an object of class
\code{gpd} are simple confidence intervals based on assumptions
of normality that are likely to be far from the truth in many cases.
A better approach would be to use profile likelihood, and we intend
to implement this method at a future date.  Alternatively, the credible intervals returned by using Bayesian estimation and the predict method for class "bgpd" will tend to give a better representation of the asymmetry of the estimated intervals around the parameter point estimates.
}
\keyword{ methods }
