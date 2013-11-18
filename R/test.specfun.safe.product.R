test.specfun.safe.product <-
function() {
  prod <- texmex:::.specfun.safe.product
  
  ## simple values
  x <- runif(10, -5, 5)
  y <- runif(length(x), -5, 5)
  
  checkEqualsNumeric(prod(x, y), pmax(x*y, -1),
                     "safe product: simple values")
  
  ## complicated values
  x <- c(0, 0, 1, 1)
  y <- c(Inf, -Inf, Inf, -Inf)
  
  res <- c(0, 0, Inf, -1)
  checkEqualsNumeric(prod(x, y), res,
                     "safe product: complicated values")
}
