context("dgev")

test_that("dgev behaves as it should", {
  evd.dgev <- .evd.dgev

  myTest <- function(mu, sig, xi, label){
      myd <- sapply(1:nreps, function(i) dgev(x[,i], mu[i], sig[i], xi[i]))
      ed <- sapply(1:nreps, function(i) evd.dgev(x[,i], loc=mu[i], scale=sig[i], shape=xi[i]))
      expect_equal(ed, myd, label=label)
  }

  set.seed(20181011)

  #*************************************************************
  # Test dgev. Note that .evd.dgev is NOT vectorized.

  nreps <- 100
  nsim <- 1000
  p <- matrix(runif(3 * nreps, -1, 1), ncol=3)
  p[, 2] <- p[, 2] + 1

  x <- sapply(1:nreps,
              function(i) rgev(nsim, mu=p[i, 1], sigma=p[i, 2], xi=p[i, 3]))

  myTest(mu=p[, 1], sig=p[, 2], xi=p[, 3], label="dgev: random xi")

  #*************************************************************
  # Test dgev when some or all of xi == 0

  p[sample(1:nreps, nreps / 2), 3] <- 0
  x <- sapply(1:nreps, function(i) rgev(nsim, mu=p[i, 1], sigma=p[i, 2], xi=p[i, 3]))
  myTest(mu=p[, 1], sig=p[, 2], xi=p[, 3], label="dgev: some zero xi")

  p[, 3] <-  0
  x <- sapply(1:nreps, function(i) rgev(nsim, mu=p[i, 1], sigma=p[i, 2], xi=p[i, 3]))
  myTest(mu=p[, 1], sig=p[, 2], xi=p[, 3], label="dgev: all zero xi")

  #*************************************************************
  # Test vectorization of dgev.

  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)
  mu <- runif(nsim, -10, 10)

  x <- rgev(nsim, mu, sig, xi)
  myd <- dgev(x, mu, sig, xi)

  ed <- sapply(1:nsim, function(i) evd.dgev(x[i], loc=mu[i], scale=sig[i], shape=xi[i]))
  expect_equal(ed, myd, label="dgev:vectorisation")

  #*************************************************************
  # 6.15 test log.d argument

  ld <- dgev(x, mu, sig, xi, log.d=TRUE)
  expect_equal(myd, exp(ld), label="dgev:logdensity")
}
)
