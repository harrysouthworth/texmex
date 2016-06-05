context("log1prel")

test_that("log1prel behaves as it should", {
  skip_on_cran()
  skip_on_travis()
  ## pull in from the namespace

  log1prel <- texmex:::.log1prel
  exprel   <- texmex:::.exprel

  ## check simple values

  x <- runif(10, -1, 5)
  log1prel.values <- log1prel(x)

  r.values <- rep.int(1, length(x))
  flag <- x!=0
  r.values[flag] <- log1p(x[flag]) / x[flag]
  expect_that(log1prel.values, equals(r.values),
                     label="log1prel: value tests")

  ## special values
  expect_that(log1prel(c(-1, 0, Inf)), equals(c(Inf,1,0)),
                     label="log1prel: special value tests")

  ## check with reference to .exprel

  x <- runif(10, -5, 5)
  y <- runif(length(x), -5, 5)

  xi <- exprel(x * y) * y

  expect_that(log1prel(xi*x)*xi, equals(y), "log1prel: exprel inversion")

}
)