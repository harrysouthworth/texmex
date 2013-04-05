qgev <- function(p, mu, sigma, xi, lower.tail=TRUE, log.p=FALSE){
  ## this code is trickier than it looks.

  ## we use G ~ (exp(-xi log E) - 1) / xi

  ## there's no argument checking in qexp, so we do this to match the
  ## existing behaviour

  if ((!log.p) && any((p < 0) || (p > 1))) {
    stop("p must lie between 0 and 1 if log.p=FALSE")
  }

  ## get the quantiles of the standard exponential
  ## change the sense of lower tail (the minus sign above)

  neg.exp.quantiles <- -log(qexp(p, lower.tail=!lower.tail,
                                 log.p=log.p))
  standard.gev <- .exprel(neg.exp.quantiles * xi) * neg.exp.quantiles

  ## and now shift and scale
  mu + sigma * standard.gev
}


test.qgev <- function() {
  ## get the probabilities that we'll use and sort them
  ## into ascending order for safekeeping
  probabilities <- sort(runif(10))

  core.sanity.test <- function(xi) {
    base.quantiles <- qgev(probabilities, 0, 1, xi)
    ## check that the values are ascending
    checkTrue(all(diff(base.quantiles) >= 0), "qgev: ascending quantiles")
    ## and check that we're descending correctly for upper tail
    bq2 <- qgev(probabilities, 0, 1, xi, lower.tail=FALSE)
    checkTrue(all(diff(bq2) <= 0), "qgev: descending quantiles")
    ## does lower.tail work
    checkEqualsNumeric(base.quantiles,
                       qgev(1 - probabilities, 0, 1, xi, lower.tail=FALSE),
                       "qgev: lower.tail works correctly")
    ## does log.p work?
    checkEqualsNumeric(base.quantiles,
                       qgev(log(probabilities), 0, 1, xi, log.p=TRUE),
                       "qgev: log.p works")
    ## check shift and scale property
    sigma <- rexp(1)
    mu    <- runif(1, -5, 5)
    shifted <- mu + sigma * base.quantiles
    checkEqualsNumeric(shifted,
                       qgev(probabilities, mu, sigma, xi),
                       "qgev: shift and scale")
  }

  lapply(c(0, seq(-5, 5, length.out=10)), core.sanity.test)

  ## known values
  checkEqualsNumeric(-log(log(2)),
                     qgev(0.5, 0, 1, 0),
                     "qgev: median match at zero xi")
  xi <- seq(-2, 2, length=10)
  checkEqualsNumeric(qgev(0.5, 0, 1, xi),
                     expm1(-log(log(2))*xi) / xi,
                     "qgev: median match at nonzero xi")
}
