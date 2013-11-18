test.gpdRangeFit <-
function(){
  par(mfrow=c(2,1))
  res <- gpdRangeFit(rain, umin=0, umax=50, nint=20)
  plot(res, pch=16, main=c("Figure 4.2 of Coles (2001)",""), addNexcesses=FALSE)
  checkEquals(names(res),c("th", "par","hi","lo","data"),msg=pst("gpdRangeFit: returned list with correct names"))
}
