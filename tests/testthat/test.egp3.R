context("egp3")

test_that("egp3 family behaves as it should", {
  library(MASS)
  rmod <- rlm(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver, method="MM", c=3.44)
  liver$r <- resid(rmod)
  
  # Check that GPD and EGP3 match when log(k) is effectively constrained to 0
  gpmod <- evm(liver[liver$dose == "D", "r"], qu=.6)
  pp <- list(c(0, 0, 0), diag(c(10^(-10), 10^4, 10^4)))
  e3mod <- evm(liver[liver$dose == "D", "r"], qu=.6, family=egp3, priorParameters=pp)

  # Test the point estimates, standard errors and t-values all at once
  expect_that(coef(summary(gpmod)), equals(coef(summary(e3mod))[-1, ], tol=.0001), label="egp3: matches gpd")

  # Check SEs on return levels - derivatives were worked out manually (by Paul),
  # also by Sage and (I think) by Wolfram
  rgp <- do.call("rbind", predict(gpmod, M=seq(100, 2000, by=100), se.fit=TRUE))
  re3 <- do.call("rbind", predict(e3mod,  M=seq(100, 2000, by=100), se.fit=TRUE))
  expect_that(rgp, equals(re3, tol=.0001), label="egp3: return levels")
})
