test_that("spgpd family does what it ought", {
  library(tidyverse)
  library(texmex)
  library(gridExtra)

  ## Simulate a large number of values so that we know we'll get close to the truth
  set.seed(20201216)
  n <- 5000
  simdat <- data.frame(i = 1:n, x = runif(n, .25, .5)) %>%
    mutate(y = rgpd(n, xi = x, sigma = 1))

  g <- evm(y, data = simdat, th = 0, xi = ~ x)
  s <- evm(y, data = simdat, th = 0, xi = ~ x, family = spgpd)

  gc <- coef(g)
  gc[length(gc)] <- log(gc[length(gc)])
  sl <- spgpd$log.lik(g$data, 0)(gc)
  gl <- gpd$log.lik(g$data, 0)(coef(g))
  expect_equal(sl, gl, label = "spgpd: log.lik returns same as gpd$log.lik")

  gc <- coef(g)
  b <- seq(gc[3] - .2, gc[3] + .2, length.out = 50)
  gl <- rep(0, 50) -> sl
  for (i in 1:50){
    gc[3] <- b[i]
    gl[i] <- gpd$log.lik(g$data, 0)(gc)
    gc[3] <- log(b[i])
    sl[i] <- spgpd$log.lik(s$data, 0)(gc)
  }

  expect_equal(gl, sl, label = "spgpd: profiling log.lik over beta same as gpd$log.lik")

  #plot(b, gl)
  #plot(b, sl)

  expect_equal(cor(resid(g), resid(s)), 1, tolerance = 1e-4, label = "spgpd: gpd and spgpd residuals are identical")

  ################################ Return levels ###############################
  n <- 500
  simdat <- data.frame(i = 1:n, x = runif(n, .25, .5)) %>%
    mutate(y = rgpd(n, xi = x, sigma = 1))

  g <- evm(y, data = simdat, th = 0, xi = ~ x)
  s <- evm(y, data = simdat, th = 0, xi = ~ x, family = spgpd)

  pg <- predict(g, M = seq(100, 1000, by = 100), ci.fit = TRUE)
  ps <- predict(s, M = seq(100, 1000, by = 100), ci.fit = TRUE)

  expect_identical(pg$obj$M.100[, "x"], ps$obj$M.100[, "x"],
                   label = "spgdp: return levels for evmOpt computed at same values as for gpd")
  expect_equal(pg$obj$M.100[, "RL"], ps$obj$M.100[, "RL"], tolerance = .01,
                   label = "spgdp: return levels for evmOpt are similar to those from gpd")


  ################################# Fit by MCMC ################################
  set.seed(20201217)
  n <- 5000
  simdat <- data.frame(i = 1:n, x = runif(n, .25, .5)) %>%
    mutate(y = rgpd(n, xi = x, sigma = 1))

  bg <- evm(y, data = simdat, th = 0, xi = ~ x, method = "sim")
  bs <- evm(y, data = simdat, th = 0, xi = ~ x, family = spgpd, method = "sim")

  ## Check posterior means. coef returns the exponent of the mean, so work with linear predictors
  bgp <- cbind(bg$param[, 1], bg$param[, 2] + simdat$x * bg$param[, 3])
  bsp <- cbind(bs$param[, 1], bs$param[, 2] + simdat$x * exp(bs$param[, 3]))

  #par(mfrow = c(2, 2))
  #hist(bgp[, 1]); hist(bgp[, 2])
  #hist(bsp[, 1]); hist(bsp[, 2])
  expect_equal(colMeans(bgp), colMeans(bsp), tolerance = .01,
               label = "spgpd: posterior means of linear predictors are similar")
  expect_equal(apply(bgp, 2, sd), apply(bsp, 2, sd), tolerance = .01,
               label = "spgpd: posterior standard deviations of linear predictors are similar")

  n <- 500
  simdat <- data.frame(i = 1:n, x = runif(n, .25, .5)) %>%
    mutate(y = rgpd(n, xi = x, sigma = 1))

  bg <- evm(y, data = simdat, th = 0, xi = ~ x, method = "sim")
  bs <- evm(y, data = simdat, th = 0, xi = ~ x, family = spgpd, method = "sim")

  bgpred <- linearPredictors(bg)
  bspred <- linearPredictors(bs)

  bgpred <- predict(bg, M = seq(500, 9000, by = 100), ci.fit = TRUE)
  bspred <- predict(bs, M = seq(500, 9000, by = 100), ci.fit = TRUE)

  par(mfrow = c(2, 3))
  plot(bgpred)
  plot(bspred)

})
