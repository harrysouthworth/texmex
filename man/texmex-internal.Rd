\name{texmex-internal}
\alias{hist.evmOpt}
\alias{qqevm}
\alias{ppevm}
\alias{qgpd2}
\alias{plotrl.evm}
\alias{plotrl.evmOpt}
\alias{u2gpd}
\alias{revTransform}
\alias{mexTransform}
\alias{evmFit}
\alias{gpd.info}
\alias{ConstraintsAreSatisfied}
\alias{PosGumb.Laplace.negloglik}
\alias{PosGumb.Laplace.negProfileLogLik}
\alias{getPlotRLdata}
\alias{plotRLevm}
\alias{rFrechet}
\alias{rMaxAR}

\alias{addCoefficients}
\alias{addCovariance}
\alias{constructEVM}

\alias{texmexPst}
\alias{namesBoot2sim}

\title{Internal functions for texmex}
\description{Internal functions used by the \code{texmex} package.}

\usage{
\method{hist}{evmOpt}(x, xlab, ylab, main, ...)
qqevm(object, nsim = 1000, alpha = 0.05)
ppevm(object, nsim = 1000, alpha = 0.05)
qgpd2(N, sigma = 1, xi = 1, u = 0, la = 1)
u2gpd(u, p=1, th=0, sigma, xi)
mexTransform(x, method = "mixture", divisor = "n+1", na.rm=TRUE,
     margins="laplace")
revTransform(x, data, qu, th=0, sigma=1, xi=0, method="mixture")
evmFit(data, family, ..., prior="none", start=NULL,
        priorParameters = NULL, maxit = 10000, trace = 0, hessian = TRUE)
gpd.info(o, method="observed")
ConstraintsAreSatisfied(a,b,z,zpos,zneg,v)
PosGumb.Laplace.negloglik(yex, ydep, a, b, m, s, constrain, v, aLow)
PosGumb.Laplace.negProfileLogLik(yex, ydep, a, b, constrain, v, aLow)
namesBoot2sim(bootobject)
getPlotRLdata(object, alpha, RetPeriodRange)
plotRLevm(M,xm,polycol,cicol,linecol,ptcol,n,xdat,pch,smooth,xlab,ylab,
  main,xrange,yrange)
plotrl.evmOpt(object, alpha = 0.05, xlab, ylab, main, pch = 1, ptcol
                 = 2, cex = 0.75, linecol = 4, cicol = 0, polycol = 15,
                 smooth = FALSE, RetPeriodRange = NULL)
rFrechet(n)
rMaxAR(n,theta)
addCoefficients(o)
addCovariance(o, family, cov)
constructEVM(o, family, th, rate, prior, modelParameters, call,
            modelData, data, priorParameters, cov)
texmexPst(msg, Family)
}

\arguments{
\item{x, object, data}{Object to be used by plot functions, vector to be converted.}
\item{xlab, ylab, main, pch, cex, linecol, polycol,
  RetPeriodRange, ...}{Arguments to plot functions.}
\item{N}{ Number of observations corresponding to N-observation
  return level to be calculated.}
  \item{modelData}{List containing response data and design matrices.}
  \item{la}{Rate at which threshold is exceeded.}
  \item{alpha}{Control nominal coverage of condfidence intervals.
    Defaults to \code{alpha = 0.05}.}
\item{nsim}{Number of simulated datasets to use in computing confidence intervals.}
\item{u}{Uniform deviates to be converted to GPD deviates.}
\item{p, th, qu, sigma, xi}{Parameters of GPD distribution.}
\item{method}{Argument to \code{mexTransform}: how to convert.
  When \code{method = "mixture"}, the upper tail of the
  distribution is modelled using a generalized Pareto distribution and the remainder
  is approximated using the empirical distribution.  Also argument
  to \code{gpd.info} which currently does nothing.}
\item{divisor}{Divisor used in estimation of empirical distribution.}
\item{na.rm}{Whether or not to remove missing values.}
\item{margins}{Form of margins to which to transform \code{x}. Can
  take values \code{margins="laplace"} or \code{margins="gumbel"}.}
\item{start, priorParameters, maxit, trace}{Arguments
    supplied to \code{gpd}, \code{migpd} or \code{mex}, or inferred from
    those functions after some preprocessing.}
\item{hessian}{Argument passed to \code{optim}. Logical.}
\item{o}{An object of class 'evmOpt', or one containing elements of such
an object.}
\item{a,b}{Dependence parameters of the Heffernan and Tawn dependence model.}
\item{m,s}{Nuisance parameters of the Heffernan and Tawn dependence model.}
\item{z,zpos,zneg}{Quantiles of the residuals under the fitted Heffernan
  and Tawn model, asymptotic positive dependence and asymptotic negative
  dependence respectively.}
\item{v}{Positive scalar, tuning parameter for constrained estimation
  of Heffernan and Tawn dependence model under estimation with Laplace
  marginal distributions.}
\item{constrain}{Logical. Whether to carry out estimation of Heffernan
  and Tawn model parameters under correct stochastic ordering of fitted
  model and asymptotic positive/negative dependence models.}
\item{aLow}{Lower bound for dependence parameter \code{a}. This depends
  on the marginal distribution under which the dependnece model is
  being fittted.  Under Gumbel margins, the lower bound is 0 and under
  Laplace margins, the lower bouind is -1.}
\item{yex, ydep}{Data for model estimation: \code{yex} is the
  explanatory variable on which the model conditions, and \code{ydep}
  is the dependent variable.}
\item{bootobject}{ Argument to \code{namesBoot2bgpd} which restructures
  an object of class evmBoot to resemble one of class evmSim, which can
  then use methods for the evmSim class.}
\item{M,xm,cicol,ptcol,n,xdat,smooth,xrange,yrange}{Arguments to
  plotRLgpd which is a worker function, does the actual plotting
  for \code{\link{plotrl.evm}}, \code{\link{plot.rl.evmSim}},
  \code{\link{plot.rl.evmBoot}}.}
\item{theta}{Argument to \code{rFrechet} and \code{rMaxAR}, the
  dependence parameter \code{theta}. Takes values between 0 and 1,
  with 0 corresponding to perfect dependence and 1 to independence.}
\item{family, rate, prior, modelParameters, call,
  cov}{Information used to obtain various components
  of an object of class 'evmOpt'.}
\item{msg, Family}{A message to print and a Family}
}

\details{
None of these functions are intended to be used explicitly.

The plotting functions are used internally by \code{\link{plot.evmOpt}}.

Some of the code is based on code that appears in the \code{ismev} package,
originally written by Stuart Coles, the evd package by Alec Stephenson
and extRemes package by Eric Gilleland, Rick Katz and Greg Young.

Code to carry out estimation of H+T2004 under Laplace margins and
constrained estimation was written by Yiannis Papastathopoulos, and
is used here for validation purposes.
}
\author{Harry Southworth, Janet E. Heffernan, Yiannis Papastathopoulos.}



