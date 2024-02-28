set.seed(20101111)
probabilities <- runif(10)

xi.values <- c(0, seq(-5, 5, length.out=10))

core.sanity.test <- function(xi) {
  randoms <- sort(rgev(10, 0, 1, xi))
  expect_true(all(diff(pgev(randoms, 0, 1, xi)) >= 0),
              info="pgev: ascending")
  expect_true(all(diff(pgev(randoms, 0,1,xi,lower.tail=FALSE)) <= 0),
              info="pgev: descending")

  expect_equal(log(pgev(randoms, 0, 1,xi)),
               pgev(randoms, 0, 1, xi, log.p=TRUE),
               info="pgev: log.p (1)")
  expect_equal(log(pgev(randoms, 0, 1, xi, lower.tail=FALSE)),
               pgev(randoms, 0, 1, xi, lower.tail=FALSE, log.p=TRUE),
               info="pgev: log.p (2)")

  mu <- runif(1, -5, 5)
  sigma <- rexp(1)

  expect_equal(pgev(randoms, 0, 1,xi),
               pgev(mu + sigma * randoms, mu, sigma, xi),
               info="pgev: shift and scale")
} # Close core.sanity.tests

qgev.comparison <- function(xi) {
  quantiles <- qgev(probabilities, 0, 1, xi)
  my.probs  <- pgev(quantiles, 0, 1, xi)
  expect_equal(my.probs, probabilities, info="pgev:straighttest")

  quantiles <- qgev(probabilities, 0, 1, xi, lower.tail=FALSE)
  my.probs  <- pgev(quantiles, 0, 1, xi, lower.tail=FALSE)
  expect_equal(my.probs, probabilities, info="pgev:lowertail")

  my.probs.2 <- pgev(quantiles, 0, 1, xi, lower.tail=TRUE)
  expect_equal(probabilities+my.probs.2,
               rep(1, length(probabilities)),
               info="pgev: tail flip")
} # Close qgev.comparison

lapply(xi.values, core.sanity.test)
lapply(xi.values, qgev.comparison)
