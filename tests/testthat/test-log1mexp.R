context("log1mexp")

test_that("simple values work", {
  skip_on_cran()
  skip_on_travis()

    x <- seq(-2, -0.01, length.out=100)
    expect_equal(.log1mexp(x), log(1-exp(x)))
    expect_equal(.log1mexp(-Inf), 0)
})

test_that("special cases work", {
  skip_on_cran()
  skip_on_travis()

    expect_true(is.na(.log1mexp(NA)))
    expect_true(is.nan(.log1mexp(0/0)))
    expect_true(is.nan(.log1mexp(Inf)))
})
