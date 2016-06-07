context("copula")

test_that("copula behaves as it should", {
  skip_on_cran()
  skip_on_travis()
    fun <- function(d) apply(d,2,function(x)(1:n)[rank(x)])/(1+n)
  n <- 200
  
  u2 <- cbind(sample(n),sample(n))
  d2 <- fun(u2)
  
  u3 <- cbind(sample(n),sample(n),sample(n))
  d3 <- fun(u3)
  
  expect_that(d2, equals(copula(u2)$copula), label="copula:2dimensional")
  expect_that(d3, equals(copula(u3)$copula), label="copula:3dimensional")
  
  op <- options()
  options(show.error.messages=FALSE)
  expect_that(copula(TRUE), throws_error(), label="copula:exception")
  expect_that(copula("text"), throws_error(), label="copula:exception")
  options(op)
}
)
