test.qgev <-
function() {
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
