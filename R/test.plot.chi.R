test.plot.chi <-
function(){
  chi <- chi(wavesurge)
  par(mfrow=c(1,2),pty="m")
  res <- plot(chi,mainChiBar="Figure 8.11 of Coles (2001)\nChi Bar")
  plot(chi, show=c("Chi"=FALSE,"ChiBar"=TRUE))
  plot(chi, show=c("Chi"=TRUE,"ChiBar"=FALSE))
  checkEquals(res,NULL,msg = "plot.chi: check successful execution")
}
