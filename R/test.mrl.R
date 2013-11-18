test.mrl <-
function(){
  par(mfrow=c(1,1))
  res <- mrl(rain)
  res <- plot(res, , main="Figure 4.1 of Coles (2001)")
  checkEquals(res, NULL, msg="mrlPlot: check execution")
}
