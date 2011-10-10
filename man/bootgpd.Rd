\name{bootgpd}
\alias{bootgpd}
\alias{print.bootgpd}
\alias{show.bootgpd}
\alias{summary.bootgpd}
\alias{summary.bootgpd}
\alias{show.summary.bootgpd}
\alias{print.summary.bootgpd}
\alias{plot.bootgpd}
\title{Parametric bootstrap for generalized Pareto models
}
\description{
Parametric bootstrap for generalized Pareto models with or without covariates
}
\usage{bootgpd(x, R = 100, trace = 10)
\method{plot}{bootgpd}(x, col="blue", border=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{An object of class \code{'gpd'}}
\item{R}{ The number of replicates to simulate. Defaults to \code{R = 100}.}
\item{trace}{Report progress to the user every \code{trace} replicates. Defaults to \code{trace = 10}.}
\item{col}{Colour of histogram bars}
\item{border}{Whether or not to draw borders between histogram bars}
\item{...}{Further arguments to \code{plot} method.}
}
\details{
The design matrices for phi and xi are held fixed and random deviates from the fitted GPD corresponding to these design matrices are generated. Especially for small sample sizes,
non-parameteric bootstrapping of GPD models can result in unreasonable
distributions (such as bimodal)
due to small numbers of observed extreme values having considerable
influence on parameter estimates in a minority of samples. Therefore,
only parametric bootstrapping is implemented here.

The \code{print} method returns the original point estimates, bootstrap means,
bootstrap estimates of bias and standard deviation. The bootstrap
median is also returned. 

The \code{summary} method returns the same as the \code{print} method, but also the
bootstrap estimate of the covariance of the parameters. When printed,
the correlation (not covariance) is displayed. The covariance might
be wanted to pass into \code{gpd} using \code{method = "simulate"}. In
particular, in some circumstances the numerical estimate of the Hessian
of the parameters can be unstable, resulting in the Metropolis algorithm
using a proposal distribution that is too narrow. This shows up as the
acceptance rate of the algorithm being too high (above about 45\%).
Then, using a bootstrap estimate might be preferable.

The \code{plot} method displays histograms and kernel density estimates.
}
\value{
  \item{call}{The function call}
  \item{replicates}{The bootstrap parameter estimates}
  \item{summary.bootgpd: margins }{Summary of the marginal parameter estiamtes}
  \item{summary.bootgpd: covariance}{Covariance of the parameter estimates}
}
\author{Harry Southworth}

\seealso{\code{\link{gpd}}}
\examples{
mod <- gpd(log(ALT.M), data=liver, qu=.7, xi=~as.numeric(dose))
bmod <- bootgpd(mod)
summary(bmod)
par(mfrow=c(1,3))
plot(bmod) 
}
\keyword{ models }

