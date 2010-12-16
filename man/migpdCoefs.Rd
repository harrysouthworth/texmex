\name{migpdCoefs}
\alias{migpdCoefs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Change values of parameters in a migpd object }
\description{
  Change the values of parameters in a \code{migpd} object. You
  might want to do this after modelling marginal distributions
  as functions of covariates.
}
\usage{
migpdCoefs(object, which, coefs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{migpd}. }
  \item{which}{ Which models in the \code{migpd} object you want to change. }
  \item{coefs}{ The coefficients that you want to change to. If \code{which}
                has length 1, \code{coefs} can be a vector of parameters.
                Otherwise, it should be a list of vectors, and the list should
                have the same length as \code{which} }
}

\value{
  A \code{migpd} object. See the help for \code{\link{migpd}}.
}

\author{ Harry Southworth }
\seealso{ \code{\link{migpd}} }
\examples{
library(MASS)
liver <- liver
liver$ndose <- as.numeric(liver$dose)
d <- data.frame(alt = resid(rlm(log(ALT_M) ~ log(ALT_B) + ndose, data=liver)),
                ast = resid(rlm(log(AST_M) ~ log(AST_B) + ndose, data=liver)),
                alp = resid(rlm(log(ALP_M) ~ log(ALP_B) + ndose, data=liver)),
                tbl = resid(rlm(log(TBL_M) ~ log(TBL_B) + ndose, data=liver)))

Dgpds <- migpd(d[liver$dose == "D", 1:4], mqu=.7)

d$ndose <- liver$ndose
galt <- gpd(alt, data=d, qu=.7, xi = ~ ndose)
gast <- gpd(ast, data=d, qu=.7, xi = ~ ndose)
galp <- gpd(alp, data=d, qu=.7, xi = ~ ndose)

altco <- c(coef(galt), coef(galt)[2] + 4 * coef(galt)[3])
astco <- c(coef(gast), coef(gast)[2] + 4 * coef(gast)[3])
alpco <- c(coef(galp), coef(galp)[2] + 4 * coef(galp)[3])

Dgpd <- migpdCoefs(Dgpds, which=c("alt", "ast", "alp"),
                   coefs=list(altco, astco, alpco))

}
\keyword{ multivariate }

