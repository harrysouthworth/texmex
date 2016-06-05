context("qgev")

test_that("qgev behaves as it should", {
  skip_on_cran()
  skip_on_travis()
function() {
  ## get the probabilities that we'll use and sort them
  ## into ascending order for safekeeping
  probabilities <- sort(runif(10))
  
  core.sanity.test <- function(xi) {
    base.quantiles <- qgev(probabilities, 0, 1, xi)
    ## check that the values are ascending
    expect_that(all(diff(base.quantiles)>=0), is_true(), "qgev:ascendingquantiles")
    ## and check that we're descending correctly for upper tail
    bq2 <- qgev(probabilities, 0, 1, xi, lower.tail=FALSE)
    expect_that(all(diff(bq2)<=0), is_true(), "qgev:descendingquantiles")
    ## does lower.tail work
    expect_that(base.quantiles, equals(), 
                qgev(1 - probabilities, 0, 1, xi, lower.tail=FALSE),
                "qgev: lower.tail works correctly")
    ## does log.p work?
    expect_that(base.quantiles, equals(), 
                qgev(log(probabilities), 0, 1, xi, log.p=TRUE),
                "qgev: log.p works")
    ## check shift and scale property
    sigma <- rexp(1)
    mu    <- runif(1, -5, 5)
    shifted <- mu + sigma * base.quantiles
    expect_that(shifted, equals(), 
                       qgev(probabilities, mu, sigma, xi),
                       "qgev: shift and scale")
  } # Close core.sanity.test
  
  lapply(c(0, seq(-5, 5, length.out=10)), core.sanity.test)
  
  ## known values
  expect_that(-log(log(2)), equals(), 
                     qgev(0.5, 0, 1, 0),
                     "qgev: median match at zero xi")
  xi <- seq(-2, 2, length=10)
  expect_that(qgev(0.5, equals(0), 1,xi),
                     expm1(-log(log(2))*xi) / xi,
                     "qgev: median match at nonzero xi")
}})
