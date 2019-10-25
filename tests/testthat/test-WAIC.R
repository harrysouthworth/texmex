context("WAIC")

testthat("WAIC, AIC and DIC are broadly in alignment on well behaved data", {
  set.seed(1234)
  b <- function(){
    i <- sample(1:nrow(liver), size=nrow(liver), replace=TRUE)
    m <- evm(log(ALT.M / ALT.B), data = liver[i, ], qu=.7, eta = ~ as.numeric(dose),
             method = "sim", iter = 2000, thin = 1, family = cgpd)
    # Use cgpd to avoid xi < -0.5 whilst bootstrapping

    AIC(m)
  }

  aic <- t(replicate(100, b()))
  expect_gte(min(cor(aic)), .98, label="WAIC: AIC, DIC and WAIC are well correlated")
})
