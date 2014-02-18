test.migpd <-
function(){
  
  # values from Heffernan and Tawn (2004) Table 4.
  # Note values in published Table 4 for u_{Xi} in cols NO2 and NO Winter were reversed.
  
  htsummer <- rbind(mqu=c(43, 43, 66.1, 22, 46),
                    mth = c(.9, .7, .7, .85, .7),
                    sigma = c(15.8, 9.1, 32.2, 42.9, 22.8),
                    xi = c(-.29, .01, .02, .08, .02))
  
  htwinter <- rbind(mqu=c(28, 49, 151.6, 23, 53),
                    mth = rep(.7, 5),
                    sigma = c(6.2, 9.3, 117.4, 19.7, 37.5),
                    xi = c(-.37, -.03, -.09, .11, -.2))
  
  summer.gpd <- summary(migpd(summer, mqu=htsummer[2,],penalty="none"),verbose=FALSE)
  winter.gpd <- summary(migpd(winter, mqu=htwinter[2,],penalty="none"),verbose=FALSE)
  
  tol <- c(1, 0.05, .5, 0.5)
  for(i in 1:4){
    checkEqualsNumeric(htsummer[i,], summer.gpd[i,],tolerance=tol[i],msg=paste("migpd: Table 4 summer",i))
    checkEqualsNumeric(htwinter[i,], winter.gpd[i,],tolerance=tol[i],msg=paste("migpd: Table 4 winter",i))
  }
  
  # check excecution for 2-d data
  
  wavesurge.fit <- migpd(wavesurge,mqu=.7)
  checkEqualsNumeric(wavesurge.fit$models$wave$loglik, evm(wavesurge$wave,qu=0.7)$loglik,
                     tolerance=0.001,msg="migpd: 2-d data gpd fit wave")
  
}
