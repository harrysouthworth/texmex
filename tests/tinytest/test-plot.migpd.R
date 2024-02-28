par(mfrow=c(2,2))
mod <- migpd(winter, mqu=.7, penalty = "none")
res <- plot(mod)
expect_equal(res, NULL, info="plot.migpd:successfulexecution")
