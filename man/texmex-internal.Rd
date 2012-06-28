\name{texmex-internal}
\alias{hist.gpd}
\alias{qqgpd}
\alias{ppgpd}
\alias{qgpd2}
\alias{plotrl.gpd}
\alias{u2gpd}
\alias{revTransform}
\alias{mexTransform}
\alias{fittedGPDscale}
\alias{fittedGPDshape}
\alias{gpdFit}
\alias{info.gpd}
\alias{ConstraintsAreSatisfied}
\alias{PosGumb.Laplace.negloglik}
\alias{PosGumb.Laplace.negProfileLogLik}
\alias{Profile_likelihood_HT_unc}
\alias{roots}
\alias{estimate_HT_KPT_joint_posneg_nm}
\alias{initial_posneg}
\alias{inv_Laplace}
\alias{Profile_likelihood_cd_nm_joint_D_KT}
\alias{Profile_likelihood_cd_nm_joint_D_KT_neg}
\alias{profile_minmax_joint_posneg_KT}
\alias{Dcond}
\alias{estimate_HT}
\alias{gpdDelta}
\alias{addCov}
\alias{namesBoot2bgpd}
\alias{plotRLgpd}
\alias{rFrechet}
\alias{rMaxAR}
\alias{test.closures}
\alias{.extRemes.decluster.intervals}
\alias{.extRemes.decluster.runs}
\alias{.extRemes.exi.intervals}
\title{Internal functions for texmex}
\description{Internal functions used by the \code{texmex} package.}

\usage{
\method{hist}{gpd}(x, xlab, ylab, main, ...)

qqgpd(object, nsim = 1000, alpha = 0.05, xlab, ylab, main, plot = TRUE,
     ylim = "auto", pch=1, col = 2, cex = 0.75, linecol = 4, intcol = 0,
     polycol = 15)
ppgpd(object, nsim = 1000, alpha = 0.05, xlab, ylab, main,
     pch=1, col = 2, cex = 0.75, linecol = 4, intcol = 0,
     polycol = 15)
qgpd2(N, sigma = 1, xi = 1, u = 0, la = 1)
u2gpd(u, p=1, th=0, sigma, xi)
mexTransform(x, method = "mixture", divisor = "n+1", na.rm=TRUE, 
     margins="laplace")
fittedGPDscale(o)
fittedGPDshape(o)
revTransform(x, data, qu, th=0, sigma=1, xi=0, method="mixture")
gpdFit(y, th, X.phi, X.xi, penalty="none", start=NULL,
        priorParameters = NULL, maxit = 10000, trace = 0, hessian = TRUE)
info.gpd(o, method="observed")
ConstraintsAreSatisfied(a,b,z,zpos,zneg,v)
PosGumb.Laplace.negloglik(yex, ydep, a, b, m, s, constrain, v, aLow)
PosGumb.Laplace.negProfileLogLik(yex, ydep, a, b, constrain, v, aLow)
Profile_likelihood_HT_unc(par,listr,x,silly=-10^(40))
estimate_HT(list,u,pars,params=TRUE)
Dcond(x,a,b,c,d,zi,zk)
roots(lev,a,c,b,d,Zj,Zk) 
Profile_likelihood_cd_nm_joint_D_KT_neg(par,listr,x,Zestfun,...,v,silly=-10^(40))
Profile_likelihood_cd_nm_joint_D_KT(par,listr,x,Zestfun,...,v,silly=-10^(40))
profile_minmax_joint_posneg_KT (pars,listdata,u,q1=0,q2=1,...,sill=-10^(40))
initial_posneg(D,...)
estimate_HT_KPT_joint_posneg_nm(pars,x,listr,params=TRUE,...,k=3)
inv_Laplace(p)
gpdDelta(A, K)
addCov(res,X)
namesBoot2bgpd(bootobject)
plotRLgpd(M,xm,polycol,cicol,linecol,ptcol,n,xdat,pch,smooth,xlab,ylab,
  main,xrange,yrange)
rFrechet(n)
rMaxAR(n,theta)
.extRemes.decluster.intervals(Z, ei) 
.extRemes.decluster.runs(Z, r) 
.extRemes.exi.intervals(Z) 
}

\arguments{
\item{x, object, data}{Object to be used by plot functions, vector to be converted.}
\item{xlab, ylab, ylim, main, pch, col, cex, linecol, intcol, polycol, ...}{Arguments to plot functions.}
  \item{N}{ Number of observations corresponding to N-observation return level to be calculated.}
  \item{la}{Rate at which threshold is exceeded.}
\item{alpha}{Control nominal coverage of condfidence intervals. Defaults to \code{alpha = 0.05}.}
\item{nsim}{Number of simulated datasets to use in computing confidence intervals.}
\item{plot}{Whether or not to plot the result.}
\item{u}{Uniform deviates to be converted to GPD deviates.}
\item{p, th, qu, sigma, xi}{Paramters of GPD distribution.}
\item{method}{Argument to \code{mexTransform}: how to convert. When \code{method = 'mixture'}, the upper tail of the
    distribution is modelled using a generalized Pareto distribution and the remainder
    is approximated using the empirical distribution.  Also argument to \code{info.gpd} which currently does nothing.}
\item{divisor}{Divisor used in estimation of empirical distribution.}
\item{na.rm}{Whether or not to remove missing values.}
\item{margins}{Form of margins to which to transform \code{x}. Can take values \code{margins="laplace"} or \code{margins="gumbel"}.}
\item{y, X.phi, X.xi, penalty, start, priorParameters, maxit, trace}{Arguments
    supplied to \code{gpd}, \code{migpd} or \code{mex}, or inferred from
    those functions after some preprocessing.}
\item{hessian}{Argument passed to \code{optim}. Logical.}
\item{o}{An
    object of class 'gpd'.}
\item{a,b}{Dependence parameters of the Heffernan and Tawn dependence model.}
\item{m,s}{Nuisance parameters of the Heffernan and Tawn dependence model.}
\item{z,zpos,zneg}{Quantiles of the residuals under the fitted Heffernan and Tawn model, asymptotic positive dependence and asymptotic negative dependence respectively.}
\item{v}{Positive scalar, tuning parameter for constrained estimation of Heffernan and Tawn dependence model under estimation with Laplace marginal distributions.}
\item{constrain}{Logical. Whether to carry out estimation of Heffernan and Tawn model parameters under correct stochastic ordering of fitted model and asymptotic positive/negative dependence models.}
\item{aLow}{Lower bound for dependence parameter \code{a}. This depends on the marginal distribution under which the dependnece model is being fittted.  Under Gumbel margins, the lower bound is 0 and under Laplace margins, the lower bouind is -1.}
\item{yex, ydep}{Data for model estimation: \code{yex} is the explanatory variable on which the model conditions, and \code{ydep} is the dependent variable.}
\item{par, listr, silly, list, pars, params, c, d, zi, zk, lev, Zj, Zk, Zestfun, 
listdata, q1, q2, sill, D, k}{Arguments to independently coded test functions implementing constrained estimation of the H+T2004 conditional dependence model with Laplace margins.}
\item{A, K}{Arguments to \code{gpdDelta} which is used internally to find
  approximate standard errors for return levels for gpd-type objects.}
\item{res,X}{ Arguments to \code{addCov} which is used internally to bind covariates on to rows of reported return levels or linear predictors.}
\item{bootobject}{ Argument to \code{namesBoot2bgpd} which restructures an object of class bootgpd to resemble one of class bgpd, which can then use methods for the bgpd class.}
\item{M,xm,cicol,ptcol,n,xdat,smooth,xrange,yrange}{Arguments to plotRLgpd which is a worker function, does the actual plotting for \code{\link{plotrl.gpd}}, \code{\link{plot.rl.gpd}}, \code{\link{plot.rl.bgpd}}.}
\item{theta}{Argument to \code{rFrechet} and \code{rMaxAR}, the dependence parameter \code{theta}. Takes values between 0 and 1, with 0 corresponding to perfect dependence and 1 to independence.}
\item{Z,r,ei}{Arguments to extRemes package functions used for testing of declustering code. Numeric Z logical indicating threshold excess, r run length for runs declustering and ei extremal index estimate.}
}

\details{
None of these functions are intended to be used explicitly.

The plotting functions are used internally by \code{\link{plot.gpd}}.

Some of the code is based on code that appears in the \code{ismev} package,
originally written by Stuart Coles, the evd package by Alec Stephenson and extRemes package by Eric Gilleland, Rick Katz and Greg Young. 

Code to carry out estimation of H+T2004 under Laplace margins and constrained estimation was written by Yiannis Papastathopoulos, and is used here for validation purposes.
}
\author{Harry Southworth, Janet E. Heffernan, Yiannis Papastathopoulos.}



