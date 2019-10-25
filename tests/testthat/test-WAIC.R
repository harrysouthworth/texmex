context("WAIC")

testthat("WAIC, AIC and DIC are broadly in alignment on well behaved data", {
  # Use priors because when xi < -0.5, I'm not sure WAIC will make any sense,
  # and maybe not even DIC
  set.seed(1234567)
  pp <- list(c(0, 0, 0), diag(c(10^4, .125, .125)))
  b <- function(){
    i <- sample(1:nrow(liver), size=nrow(liver), replace=TRUE)
    m <- evm(log(ALT.M / ALT.B), data = liver[i, ], qu=.7, xi = ~ as.numeric(dose),
             method = "sim", iter = 2000, thin = 1, priorParameters = pp)

    s <- sum((m$param[, 3] - mean(m$param[, 3])^3))

    c(AIC(m), s)
  }

  aic <- t(replicate(100, b()))
  expect_gte(min(cor(aic[, 1:3])), .98, label="WAIC: AIC, DIC and WAIC are well correlated")
})
