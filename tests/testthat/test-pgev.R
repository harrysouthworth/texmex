context("pgev")

test_that("pgev behaves as it should", {
  set.seed(20101111)
  probabilities <- runif(10)

  xi.values <- c(0, seq(-5, 5, length.out=10))

  core.sanity.test <- function(xi) {
    randoms <- sort(rgev(10, 0, 1, xi))
    expect_true(all(diff(pgev(randoms, 0, 1, xi)) >= 0),
                label="pgev: ascending")
    expect_true(all(diff(pgev(randoms, 0,1,xi,lower.tail=FALSE)) <= 0),
                label="pgev: descending")

    expect_equal(log(pgev(randoms, 0, 1,xi)),
                 pgev(randoms, 0, 1, xi, log.p=TRUE),
                 label="pgev: log.p (1)")
    expect_equal(log(pgev(randoms, 0, 1, xi, lower.tail=FALSE)),
                 pgev(randoms, 0, 1, xi, lower.tail=FALSE, log.p=TRUE),
                 label="pgev: log.p (2)")

    mu <- runif(1, -5, 5)
    sigma <- rexp(1)

    expect_equal(pgev(randoms, 0, 1,xi),
                 pgev(mu + sigma * randoms, mu, sigma, xi),
                 label="pgev: shift and scale")
  } # Close core.sanity.tests

  qgev.comparison <- function(xi) {
    quantiles <- qgev(probabilities, 0, 1, xi)
    my.probs  <- pgev(quantiles, 0, 1, xi)
    expect_equal(my.probs, probabilities, label="pgev:straighttest")

    quantiles <- qgev(probabilities, 0, 1, xi, lower.tail=FALSE)
    my.probs  <- pgev(quantiles, 0, 1, xi, lower.tail=FALSE)
    expect_equal(my.probs, probabilities, label="pgev:lowertail")

    my.probs.2 <- pgev(quantiles, 0, 1, xi, lower.tail=TRUE)
    expect_equal(probabilities+my.probs.2,
                 rep(1, length(probabilities)),
                 label="pgev: tail flip")
  } # Close qgev.comparison

  lapply(xi.values, core.sanity.test)
  lapply(xi.values, qgev.comparison)

  ##############################################################################
  # The remaining tests copy evd::pgev and run tests similar to dgev (which was
  # copied form pgev)

  evd.pgev <- .evd.pgev

  myTest <- function(mu, sig, xi, msg){
    myp <- sapply(1:nreps, function(i) pgev(x[,i], mu[i], sig[i], xi[i]))
    ep <- sapply(1:nreps, function(i) evd.pgev(x[,i], loc=mu[i], scale=sig[i], shape=xi[i]))
    expect_equal(ep, myp, label=msg)
  } # Close myTest

  set.seed(20181011)

  #*************************************************************
  # Test pgev. Note that .evd.pgev is NOT vectorized.

  nreps <- 100
  nsim <- 1000
  p <- matrix(runif(3 * nreps, -1, 1), ncol=3)
  p[, 2] <- p[, 2] + 1

  x <- sapply(1:nreps, function(i) rgev(nsim, mu=p[i, 1], sigma=p[i, 2], xi=p[i, 3]))

  myTest(mu=p[, 1], sig=p[, 2], xi=p[, 3], msg="pgev: random xi")

  #*************************************************************
  # Test pgev when some or all of xi == 0

  p[sample(1:nreps, nreps / 2), 3] <- 0
  x <- sapply(1:nreps, function(i)
    rgev(nsim, mu=p[i, 1], sigma=p[i, 2], xi=p[i, 3]))
  myTest(mu=p[, 1], sig=p[, 2], xi=p[, 3], msg="pgev: some zero xi")

  p[, 3] <-  0
  x <- sapply(1:nreps,function(i) rgev(nsim, mu=p[i, 1], sigma=p[i, 2], xi=p[i, 3]))
  myTest(mu=p[, 1], sig=p[, 2], xi=p[, 3], msg="pgev: all zero xi")

  #*************************************************************
  # Test vectorization of pgev.

  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)
  mu <- runif(nsim, -10, 10)

  x <- rgev(nsim, mu, sig, xi)
  myp <- pgev(x, mu, sig, xi)

  ep <- sapply(1:nsim, function(i)
    evd.pgev(x[i], loc=mu[i], scale=sig[i], shape=xi[i]))
  expect_equal(ep, myp, label="pgev:vectorisation")

  #*************************************************************
  # Test log.p argument

  lp <- pgev(x, mu, sig, xi, log.p=TRUE)
  expect_equal(myp, exp(lp), label="pgev:logprobabilities")

  #*************************************************************
  # Test lower tail argument

  sp <- pgev(x, mu, sig, xi, lower.tail=FALSE)
  expect_equal(myp, 1-sp, label="pgev:lowertail")

  ## check pgev when xi < 0 and value above upper limit

  xi <- -2.3
  upperProb <- pgev(-2/xi, mu=0, sigma=1, xi, lower.tail=TRUE)
  expect_equal(upperProb, 1, label="pgev:negativexi(1)")
}
)
