\name{bootgpd}
\alias{bootgpd}
\alias{print.bootgpd}
\alias{show.bootgpd}
\alias{summary.bootgpd}
\alias{summary.bootgpd}
\alias{plot.bootgpd}
\title{Parametric bootstrap for generalized Pareto models
}
\description{
Parametric bootstrap for generalized Pareto models with or without covariates
}
\usage{bootgpd(x, R = 99, trace = 10)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class 'gpd'
}
  \item{R}{ The number of replicates to run. Defaults to \code{R = 99}.
}
  \item{trace}{
Report progress to the user every \code{trace} replicates. Defaults to \code{trace = 10}.}
\item{x){An object of class \code{gpd}}
\item{object){An object of class \code{gpd}}
\item{col}{Colour of histogram bars}
\item{border}{Whether or not to draw borders between histogram bars}

}
\details{
The design matrices for phi and xi are held intact and random numbers are
generated as the bootstrap response variable. Especially for small sample sizes,
non-parameteric bootstrapping of GPD models can result in unreasonable distributions
due to small numbers of observed extreme values having considerable influence
on parameter estimates. Therefore, only parametric bootstrapping is implemented
here.
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

