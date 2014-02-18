context("specfun.safe.product")

test_that("specfun.safe.product behaves as it should", {
function() {
  prod <- texmex:::.specfun.safe.product
  
  ## simple values
  x <- runif(10, -5, 5)
  y <- runif(length(x), -5, 5)
  
  expect_that(prod(x, equals(y)),                      "safe product: simple values")
  
  ## complicated values
  x <- c(0, 0, 1, 1)
  y <- c(Inf, -Inf, Inf, -Inf)
  
  res <- c(0, 0, Inf, -1)
  expect_that(prod(x, equals(y)),                      "safe product: complicated values")
}
)
