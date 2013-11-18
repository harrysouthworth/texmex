\name{internal}
\alias{.extRemes.decluster.intervals}
\alias{.extRemes.decluster.runs}
\alias{.extRemes.exi.intervals}
\alias{Profile_likelihood_HT_unc}
\alias{roots}
\alias{Profile_likelihood_cd_nm_joint_D_KT}
\alias{Profile_likelihood_cd_nm_joint_D_KT_neg}
\alias{profile_minmax_joint_posneg_KT}
\alias{Dcond}
\alias{estimate_HT_KPT_joint_posneg_nm}
\alias{initial_posneg}
\alias{inv_Laplace}
\alias{estimate_HT}

\title{Internal functions used to run tests on texmex}

\usage{
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
  .extRemes.decluster.intervals(Z, ei)
  .extRemes.decluster.runs(Z, r)
  .extRemes.exi.intervals(Z)
}

\arguments{
  \item{par, listr, silly, list, pars, params, c, d, zi, zk, lev}{Arguments to
        independently coded test functions implementing constrained estimation
        of the H+T2004 conditional dependence model with Laplace margins.}
  \item{Zj, Zk, Zestfun, listdata, q1, q2, sill, D, k}{Arguments to
        independently coded test functions implementing constrained estimation
        of the H+T2004 conditional dependence model with Laplace margins.}
  \item{Z,r,ei}{Arguments to extRemes package functions used for testing
  of declustering code. Numeric Z logical indicating threshold excess, r
  run length for runs declustering and ei extremal index estimate.}
}

\author{
  Janet E. Heffernan and Harry Southworth
}
