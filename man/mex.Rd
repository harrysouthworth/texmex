\name{mex}
\alias{mex}
\alias{plot.mex}
\alias{predict.mex}
\alias{summary.predict.mex}
\alias{plot.predict.mex}
\title{Conditional multivariate extreme values modelling}
\description{Fit the conditional multivariate extreme value model of Heffernan and Tawn}
\usage{
mex(data, which, mth, mqu, dqu, margins = "laplace", constrain = TRUE, v = 10, 
    penalty = "gaussian", maxit = 10000, trace = 0, verbose = FALSE, priorParameters = NULL)

\method{plot}{mex}(x, quantiles = seq(0.1, by = 0.2, len = 5), col = "grey", ...)
\method{predict}{mex}(object, which, pqu=0.99, nsim=1000, trace=10, ...)
\method{summary}{predict.mex}(object, mth, probs=c(0.05, 0.5, 0.95), ...)
\method{plot}{predict.mex}(x, pch=c(1, 3), col=c(2, 8), cex=c(1, 1), ask=TRUE, ...)
}

\arguments{
  \item{data}{A numeric matrix or data.frame, the columns of which are to be modelled.}
  \item{which}{The variable on which to condition.  This can be either scalar, indicating the column number of the conditioning variable, or character, giving the column name of the conditioning variable.}
  \item{mth}{In \code{mex}, the threshold above which to fit generalized Pareto distributions.
            If this is a vector of length 1, the same threshold will be used for
            each variable. Otherwise, it should be a vector whose length is
            equal to the number of columns in \code{data}.

            In \code{summary.predict.mex}, the thresholds over which to simulate
            data from the fitted multivariate model. If not supplied, it is
            taken to be the thresholds that were used to fit the dependence model on the
            scale of the original data.}
  \item{mqu}{As an alternative to specifying the marginal GPD fitting thresholds via mth, you can
            specify the quantile (a probability) above which to fit generalized Pareto
            distributions. If this is a vector of length 1, the same quantile
            will be used for each variable. Otherwise, it should be a vector
            whose length is equal to the number of columns in data.}
  \item{dqu}{Used to specify the quantile at which to
            threshold the conditioning variable data when estimating the dependence
            parameters.  For example \code{dqu=0.7} will result in the data with the highest 30% of values of the conditioning variable being used to estimated the dependence parameters. The same
            threshold will be used for each dependent variable.  If not supplied then the default is to set \code{dqu=mqu[which]} the quantile corersponding to the threshold used to fit the marginal model to the tail of the conditioning variable.  Note that there is no requirement for the quantiles used for marginal fitting (\code{mqu}) and dependence fitting (\code{dqu}) to be the same, or for them to be ordered in any way.}
  \item{margins}{See documentation for \code{\link{mexDependence}}.}
  \item{constrain}{See documentation for \code{\link{mexDependence}}.}
  \item{v}{See documentation for \code{\link{mexDependence}}.}
  \item{penalty}{How to penalize the likelihood when estimating the marginal
            generalized Pareto distributions. Defaults to ``gaussian''. See the
            help file for \code{\link{gpd}} for more information.}
  \item{maxit}{The maximum number of iterations to be used by the optimizer.
            defaults to \code{maxit = 10000}.}
  \item{trace}{Passed internall to \code{\link{optim}}. Whether or not to inform the user of the progress of the
            optimizer. Defaults to 0, indicating no trace.}
  \item{verbose}{Whether or not to keep the user informed of progress. Defaults
            to \code{verbose = FALSE}.}
  \item{priorParameters}{Parameters of prior/penalty used for estimation of the 
            GPD parameters.  This is only used if \code{penalty = 'gaussian'}.  
            It is a named list, each element of which contains two
            components: the first component should be a vector of length 2 corresponding
            to the location of the Gaussian distribution; the second a 2x2 matrix corresponding to the covariance matrix of the distribution.
            The names should match the names of the columns of \code{data}. 
            If not provided, the default priors are independent normal, centred at zero, with 
            variance 10000 for log(sigma) and 0.25 for xi. See the details section.}
  \item{quantiles}{A vector of quantiles taking values between 0 and 1 specifying 
            the quantiles of the conditional distributions which will be plotted.}
  \item{col}{In \code{plot} method for objects of class \code{mex}, the color for points on scatterplots of residuals and original data respectively.
             In \code{plot} method for objects of class \code{predict.mex}, the colours of points for observed and
             simulated data respectively.}
  \item{x, object}{Object of class \code{mex} as returned by function \code{mex}.}
  \item{pqu}{Argument to \code{predict} method. The quantile of the conditioning 
            variable above which it will be simulated for importance sampling based prediction. Defaults to \code{pqu = .99 }. }
  \item{nsim}{Argument to \code{predict} method. The number of simulated observations 
            to be generated for prediction. }
  \item{probs}{In \code{summary} method for objects of class \code{predict.mex}: the quantiles of the conditional 
            distribution(s) to calculate.  Defaults to 5\%, 50\% and 95\%.}
  \item{pch, cex}{Plotting characters: colours and symbol expansion. The
            observed and simulated data are plotted using different symbols, controlled
            by these arguments and \code{col}, each of which should be of length 2.}
  \item{ask}{Whether or not to ask before changing the plot. Defaults to
            \code{ask = TRUE}.}
  \item{...}{Further arguments to be passed to methods.}
}

\details{The function \code{mex} works as follows. First, Generalized Pareto distributions (GPD) are fitted to the upper tails of each of the marginal distributions of the data: the GPD parameters are estimated for
         each column of the data in turn, independently of all other columns.
         Then, the conditional multivariate approach of Heffernan and Tawn is
         used to model the dependence between variables. The returned object is of class 'mex'.
         
         This function is a wrapper for calls to \code{\link{migpd}} and \code{\link{mexDependence}}, which estimate parameters of the marginal and dependence components of the Heffernan and Tawn model respectively.  See documentation of these functions for details of modelling issues including the use of penalties / priors, threshold choice and checking for convergence of parameter estimates.

     The \code{plot} method produces diagnostic plots for the fitted dependence model described by Heffernan and Tawn, 2004.  The plots are best viewed by using the plotting area split by \code{par(mfcol=c(.,.))} rather than \code{mfrow}, see examples below.  Three diagnostic plots are produced for each dependent variable:
    
    1) Scatterplots of the residuals Z from the fitted model of Heffernan and Tawn (2004) are 
    plotted against the quantile of the conditioning variable, with a lowess curve showing the local 
    mean of these points.  2) The absolute value of \code{Z-mean(Z)} is also plotted, again with the lowess curve showing 
    the local mean of these points.  Any trend in the location or scatter of these variables with the conditioning variable
    indicates a violation of the model assumption that the residuals Z are indepenendent of the conditioning
    variable.  This can be indicative of the dependence threshold used being too low. 3) The final plots show the original data (on the original scale) and the fitted quantiles (specified by \code{quantiles}) of the conditional distribution 
    of each dependent variable given the conditioning variable.  A model that fits well will have good agreement between the
    distribution of the raw data (shown by the scatter plot) and the fitted quantiles. Note that the raw data are a sample from the joint distribution, whereas the quantiles are those of the estimated conditional distribution given the value of the conditioning variable, and while these two distributions should move into the same part of the sample space as the conditioning variable becomes more extreme, they are not the same thing!

    The \code{predict} method for \code{mex} works as follows. The returned object 
         has class 'predict.mex'. Simulated values of the dependent variables are created, 
         given that the conditioning variable is above its 100\code{pqu} quantile.  
         If \code{predict.mex} is passed an object \code{object} of class \code{"mex"} 
         then the simulated values are based only on the point estimate of the dependence 
         model parameters, and the original data.  If \code{predict.mex} is passed an object 
         \code{object} of class \code{"bootmex"} then the returned value additionally 
         contains simulated replicate data sets corresponding to the bootstrap model parameter 
         estimates.  In both cases, the simulated values based on the original data
         and point estimates appear in component \code{object$data$simulated}. The
         simulated data from the bootstrap estimates appear in \code{object$replicates}.
         
         The \code{plot} method for class \code{"oredict.mex"} displays both the original 
         data and the simulated data generated above the threshold for prediction; it 
         shows the threshold for prediction (vertical line) and also the curve joining equal quantiles 
         of the marginal distributions -- this is for reference: variables that 
         are perfectly dependent will lie exactly on this curve.
         
         

}

\value{
A call to \code{mex} returns an list of class \code{mex} containing the following three items:
  \item{margins}{An object of class \code{\link{migpd}}.}
  \item{dependence}{An object of class \code{\link{mexDependence}}.}
  \item{call}{This matches the original function call.}
There are \code{plot}, \code{summary}, \code{coef} and \code{predict} methods for this class.

A call to \code{predict.mex} does the importance sampling for prediction, and returns a list of class \code{"predict.mex"} for which there are print and plot methods available.  The summary method for this class of object is intended to be used following a call to the predict method, to estimate quantiles or probabilities of threshold excesses for the fitted conditional distributions given the conditioning variable above the threshold for prediction.  See examples below.

There are \code{print}, \code{summary} and
         \code{plot} methods available for the class 'predict.mex'. 

}
\references{J. E. Heffernan and J. A. Tawn, A conditional approach for
            multivariate extreme values, Journal of the Royal Statistical
            Society B, 66, 497 - 546, 2004}
\author{
Harry Southworth, Janet E. Heffernan
}
\note{The package \code{texmex} is equipped to fit GPD models to the upper marginal tails only, not the lower tails.  This is appropriate for extrapolating into the tails of any dependent variable when dependence between this variable and the conditioning variable is positive.  In the case of negative dependence between the conditioning variable and any dependent variable, estimation of the conditional distribution of the dependent variable for extreme values of the conditioning variable would naturally visit the lower tail of the dependent variable.  Extrapolation beyond the range of the observed lower tail is not supported in the current version of \code{texmex}. In cases where negative dependence is observed and extrapolation is required into the lower tail of the dependent variable, the situation is trivially resolved by working with a reflection of the dependent variable (Y becomes -Y and so the upper and lower tails are swapped). Results can be calculated for the reflected variable then reflected back to the correct scale.   This is satisfactory when only the pair of variables (the conditioning and single dependent variable) are of interest, but when genuine multivariate (as opposed to simply bivariate) structure is of interest, this approach will destroy the dependence structure between the reflected dependent variable and the remaining dependent variables.}
\seealso{\code{\link{migpd}}, \code{\link{mexDependence}}, \code{\link{bootmex}}}
\examples{
w <- mex(winter, mqu=.7)
w
par(mfcol=c(3, 2))
plot(w)

par(mfcol=c(2,2))
p <- predict(w)
summary(p)
summary(p,probs=c(0.01,0.2,0.5,0.8,0.99))
summary(p,probs=0.5,mth=c(40,50,150,25,50))
p
plot(p)

}
\keyword{ models }
\keyword{ multivariate }
