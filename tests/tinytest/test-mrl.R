par(mfrow=c(1,1))
res <- mrl(rain)
res <- plot(res, main="Figure 4.1 of Coles (2001)")
expect_equal(res, NULL, label="mrlPlot:checkexecution")
