context("plot.migpd")

test_that("plot.migpd behaves as it should", {
    par(mfrow=c(2,2))
  mod <- migpd(winter, mqu=.7, penalty = "none")
  res <- plot(mod)
  expect_that(res, equals(NULL), }
)
