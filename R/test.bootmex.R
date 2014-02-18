test.bootmex <-
function(){ # this is a weak test - it tests the structure
  # of the output but not the correctness of the bootstrap coefficients; it will
  # also catch ERRORs (as opposed to FAILUREs) if the code breaks.  For strong
  # testing of this function, run test.predict.mex
  
  set.seed(20120327)
  smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
  wmarmod <- migpd(winter, mqu=.7,  penalty="none")
  
  mySdep <- mexDependence(smarmod, dqu=.7, which=1)
  myWdep <- mexDependence(wmarmod, dqu=.7, which=1)
  
  R <- 20
  
  mySboot <- bootmex(mySdep, R=R,trace=R+1)
  myWboot <- bootmex(myWdep, R=R,trace=R+1)
  
  checkEqualsNumeric(coef(mySdep)[[2]], mySboot$simpleDep, msg="bootmex: summer simpleDep from call with mar model")
  checkEqualsNumeric(coef(myWdep)[[2]], myWboot$simpleDep, msg="bootmex: winter simpleDep from call with mar model")
  
  checkEqualsNumeric(coef(smarmod), coef(mySboot$simpleMar), msg="bootmex: summer simpleMar from call with mar model")
  checkEqualsNumeric(coef(wmarmod), coef(myWboot$simpleMar), msg="bootmex: winter simpleMar from call with mar model")
  
  checkEqualsNumeric(dim(mySboot$boot[[1]]$Z)[2], dim(mySdep$dependence$Z)[2], msg="bootmex: summer dim of residuals from call with mar model")
  checkEqualsNumeric(dim(myWboot$boot[[1]]$Z)[2], dim(myWdep$dependence$Z)[2], msg="bootmex: winter dim of residuals from call with mar model")
  
  checkEqualsNumeric(dim(mySboot$boot[[1]]$dependence), dim(coef(mySdep)[[2]]), msg="bootmex: summer dim of coefficients from call with mar model")
  checkEqualsNumeric(dim(myWboot$boot[[1]]$dependence), dim(coef(myWdep)[[2]]), msg="bootmex: winter dim of coefficients from call with mar model")
  
  checkEqualsNumeric(R, length(mySboot$boot),msg="bootmex: number of bootstrap samples")
  checkEqualsNumeric(R, length(myWboot$boot),msg="bootmex: number of bootstrap samples")
  
  checkEqualsNumeric(dim(summer),dim(mySboot$boot[[1]]$Y),msg="bootmex: size of bootstrap data set")
  checkEqualsNumeric(dim(winter),dim(myWboot$boot[[5]]$Y),msg="bootmex: size of bootstrap data set")
  
  smexmod <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="gumbel",constrain=FALSE, which=1)
  wmexmod <- mex(winter, mqu=.7,  penalty="none", margins="gumbel",constrain=FALSE, which=1)
  
  mySboot <- bootmex(smexmod, R=R,trace=R+1)
  myWboot <- bootmex(wmexmod, R=R,trace=R+1)
  
  checkEqualsNumeric(coef(smexmod)[[2]], mySboot$simpleDep, msg="bootmex: summer simpleDep from call with mex model")
  checkEqualsNumeric(coef(smexmod)[[1]], coef(mySboot$simpleMar), msg="bootmex: summer simpleMar from call with mex model")
  
  checkEqualsNumeric(coef(wmexmod)[[2]], myWboot$simpleDep, msg="bootmex: winter simpleDep from call with mex model")
  checkEqualsNumeric(coef(wmexmod)[[1]], coef(myWboot$simpleMar), msg="bootmex: winter simpleMar from call with mex model")
  
  checkEqualsNumeric(R, length(mySboot$boot),msg="bootmex: number of bootstrap samples")
  checkEqualsNumeric(R, length(myWboot$boot),msg="bootmex: number of bootstrap samples")
  
  checkEqualsNumeric(dim(summer),dim(mySboot$boot[[1]]$Y),msg="bootmex: size of bootstrap data set")
  checkEqualsNumeric(dim(winter),dim(myWboot$boot[[5]]$Y),msg="bootmex: size of bootstrap data set")
  
  smexmod.1 <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="laplace",constrain=FALSE, which=1)
  smexmod.2 <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="laplace",constrain=TRUE,v=2, which=1)
  mySboot.1 <- bootmex(smexmod.1,R=R,trace=R+1)
  mySboot.2 <- bootmex(smexmod.2,R=R,trace=R+1)
  
  checkEqualsNumeric(coef(smexmod.1)[[2]], mySboot.1$simpleDep, msg="bootmex: summer simpleDep from call with mex model, constrain=FASLE")
  checkEqualsNumeric(coef(smexmod.1)[[1]], coef(mySboot.1$simpleMar), msg="bootmex: summer simpleMar from call with mex model, constrain=FASLE")
  checkEqualsNumeric(coef(smexmod.2)[[2]], mySboot.2$simpleDep, msg="bootmex: summer simpleDep from call with mex model, v=2")
  checkEqualsNumeric(coef(smexmod.2)[[1]], coef(mySboot.2$simpleMar), msg="bootmex: summer simpleMar from call with mex model, v=2")
  
  # check execution of for 2-d data
  
  wavesurge.fit <- migpd(wavesurge,mqu=.7)
  wavesurge.mex <- mexDependence(wavesurge.fit, dqu=0.8,which=1)
  R <- 20
  
  wavesurge.boot <- bootmex(wavesurge.mex,R=R,trace=R+1)
  
  checkEqualsNumeric(dim(wavesurge.boot$boot[[1]]$Z)[2],1,msg="bootmex: execution for 2-d data")
  checkEquals(dimnames(wavesurge.boot$boot[[1]]$Z)[[2]],names(wavesurge)[2],msg="bootmex: execution for 2-d data")
  checkEqualsNumeric(length(wavesurge.boot$boot),R,msg="bootmex: execution for 2-d data")
}
