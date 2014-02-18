test.plot.migpd <-
function(){
  par(mfrow=c(2,2))
  mod <- migpd(winter, mqu=.7, penalty = "none")
  res <- plot(mod)
  checkEquals(res,NULL,msg="plot.migpd: successful execution")
}
