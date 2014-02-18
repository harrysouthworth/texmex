context("exprel")

test_that("exprel behaves as it should", {
function(x) {
  ## pull in from the namespace because RUnit doesn't give access to
  ## internal functions.
  exprel <- texmex:::.exprel
  
  ## first check some simple values
  values <- c(-Inf, NA, seq(-5, 5, length.out=10))
  r.values <- expm1(values) / values
  exprel.values <- exprel(values)
  expect_that(exprel.values, equals(r.values),                      label="exprel: value tests")
  
  
  ## the naive R code doesn't compute these right so force the
  ## computation.
  values <- c(Inf, 0)
  exprel.values <- exprel(values)
  desired.values <- c(Inf, 1)
  expect_that(exprel.values, equals(desired.values),                      label="exprel: special value tests")
  
}
)
