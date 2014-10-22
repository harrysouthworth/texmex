context("simulate")

test_that("simulate.evm functions behave as expected", {
  ##############################################################################
  ##################################### GPD ####################################
  # Compare simulate.evmOpt with rgpd
  m <- evm(rain, qu=.95)
  s1 <- simulate(m, nsim=100, seed=1234)
  set.seed(1234)
  s2 <- rgpd(100, exp(coef(m)[1]), coef(m)[2])
  expect_that(cor(s1, s2), equals(1), label="Compare simulate.evmOpt with rgpd")

  # Compare simulate.evmSim with rgpd
  m <- evm(rain, qu=.95, method="sim", iter=1500, burn=500, thin=1)
  s1 <- simulate(m, seed=1234)
  set.seed(1234)
  s2 <- rgpd(1000, exp(m$param[, 1]), m$param[, 2])
  expect_that(cor(s1, s2), equals(1), label="Compare simulate.evmSim with rgpd")
  
  # Compare simulate.evmBoot with rgpd
  m <- evm(rain, qu=.95, method="boot", R=1000)
  s1 <- simulate(m, seed=1234)
  set.seed(1234)
  s2 <- rgpd(1000, exp(m$replicates[, 1]), m$replicates[, 2])
  expect_that(cor(s1, s2), equals(1), label="Compare simulate.evmBoot with rgpd")

  ##############################################################################
  ##################################### GEV ####################################
  # Compare simulate.evmOpt with rgev
  m <- evm(SeaLevel, portpirie, family=gev)
  s1 <- simulate(m, nsim=100, seed=1234)
  set.seed(1234)
  s2 <- rgev(100, coef(m)[1], exp(coef(m)[2]), coef(m)[3])
  expect_that(cor(s1, s2), equals(1), label="Compare simulate.evmOpt with rgev")

  # Compare simulate.evmSim with rgev
  m <- evm(SeaLevel, portpirie, family=gev, method="sim", iter=1500, burn=500, thin=1)
  s1 <- simulate(m, seed=1234)
  set.seed(1234)
  s2 <- rgev(1000, m$param[, 1], exp(m$param[, 2]), m$param[, 3])
  expect_that(cor(s1, s2), equals(1), label="Compare simulate.evmSim with rgev")
  
  # Compare simulate.evmBoot with rgev
  m <- evm(SeaLevel, portpirie, family=gev, method="boot", R=1000)
  s1 <- simulate(m, seed=1234)
  set.seed(1234)
  s2 <- rgev(1000, m$replicates[, 1], exp(m$replicates[, 2]), m$replicates[, 3])
  expect_that(cor(s1, s2), equals(1), label="Compare simulate.evmBoot with rgev")
})