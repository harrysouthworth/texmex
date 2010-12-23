\name{texmex-package}
\alias{texmex-package}
\alias{texmex}
\docType{package}
\title{
Conditional multivariate extreme values modelling.
}
\description{
Conditional multivariate extreme values modelling following the approach
of Heffernan and Tawn, 2004.
}
\details{
\tabular{ll}{
Package: \tab mex\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-12-06\cr
License: \tab GPL (>=2) | BSD\cr
}
\code{\link{gpd}}: Fit generalized Pareto distributions to data, possibly with
            covariates. Use maximum likelihood estimation, maximum
            penalized likelihood estimation or simulate from the posterior
            distribution.

\code{\link{mex}}: Fit multiple, independent generalized Pareto models
            to the columns of a data set, and estimate the conditional
            dependence structure using the method of Heffernan and
            Tawn.

\code{\link{bootmex}}: Bootstrap estimation for parameters in generalized Pareto models
                and in the dependence structure.

\code{\link{validate.texmex}}: Use the \code{svUnit} package to run lots of tests to
            establish a high degree of confidence that the functions do what
            they are supposed to.
}
\author{
Harry Southworth, Janet E. Heffernan

Maintainer: Harry Southworth <harry.southworth@astrazeneca.com>

}
\references{
J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004
}
\keyword{models}
\keyword{multivariate}
\keyword{package}
\examples{
# Analyse the winter data used by Heffernan and Tawn
mymex <- mex(winter, mqu = .7, penalty="none")
plot(mymex)
myboot <- bootmex(mymex, dqu=.7, which = "NO")
plot(myboot)
mypred <- predict(myboot,  pqu=.95)
summary(mypred , probs = c( .025, .5, .975 ))

# Analyse the liver data included in the package
library(MASS) # For the rlm function

liver <- liver[liver$ALP.M > 1,] # Get rid of outlier
liver$ndose <- as.numeric(liver$dose)

alt <- resid(rlm(log(ALT.M) ~ log(ALT.B) + ndose, data=liver, method="MM"))
ast <- resid(rlm(log(AST.M) ~ log(AST.B) + ndose, data=liver, method="MM"))
alp <- resid(rlm(log(ALP.M) ~ log(ALP.B) + ndose, data=liver, method="MM"))
tbl <- resid(rlm(log(TBL.M) ~ log(TBL.B) + ndose, data=liver, method="MM"))

r <- data.frame(alt=alt, ast=ast, alp=alp, tbl=tbl)

Amex <- mex(r[liver$dose == "A",], mqu=.7)
Bmex <- mex(r[liver$dose == "B",], mqu=.7)
Cmex <- mex(r[liver$dose == "C",], mqu=.7)
Dmex <- mex(r[liver$dose == "D",], mqu=.7)

par(mfcol=c(2,3))
plot(Amex)

plot(Dmex, col="blue")

## Take a closer look at ALT

r$ndose <- liver$ndose

altmod1 <- gpd(alt, qu=.7, phi = ~ ndose, xi = ~ ndose, data=r)
altmod2 <- gpd(alt, qu=.7, phi = ~ ndose, data=r)
altmod3 <- gpd(alt, qu=.7, xi = ~ ndose, data=r)
altmod4 <- gpd(alt, qu=.7, data=r)

# Prefer model 3, with term for xi on basis of AIC

# Lines commented out to keep CRAN robots happy
#balt3 <- gpd(alt, qu=.7, xi = ~ ndose, data=r, method="simulate")
#par(mfrow=c(3,3))
#plot(balt3)

# use longer burn-in and also thin the output

#balt3 <- thinAndBurn(balt3,burn=1000,thin=5)
#plot(balt3)

# Get some simulated values for dose D

#Dphi <- balt3$param[, 1]
#Dxi <- balt3$param[, 2] + 4 * balt3$param[, 3]

#simD <- rgpd(nrow(balt3$param), sigma=exp(Dphi), xi=Dxi, u=quantile(alt, .7))

# Those are simulated residuals. Get some baselines and transform all
# to raw scale

#b <- sample(log(liver$ALT.M), size=nrow(balt3$param), replace=TRUE)
#res <- exp(b + simD)

# estimate quantiles on raw scale
#quantile(res, prob=c(.5, .75, .9, .95, .99))

# estimate proportion exceeding 3*upper limit of normal
#mean(res > 36 * 3) # 36 is the upper limit of normal for ALT
}

