# functionality all tested already in test.predict.evm, so just check output of correct format.

n <- 1000
nx <- 9
R <- 10
nm <- 20

for(Family in list(gpd,gev)){

  pst <- function(msg) texmex:::texmexPst(msg,Family=Family)
  set.seed(20130513)

  X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
  param <- switch(Family$name,GPD=X,GEV=cbind(5,X))
  th <- switch(Family$name,GPD=0,GEV=-Inf)
  X$Y <- Family$rng(n,param,list(threshold=th))
  fit <- evm(Y,data=X,phi=~a,xi=~b,th=th,family=Family)

  o <- options(warn=-1)
  boot <- evmBoot(fit,R=R,trace=100)
  options(o)

  newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
  from <- 10; to <- 500
  M <- seq(from,to,length=nm)
  pred <- predict(boot,newdata=newX,M=M,ci=TRUE)

  expect_equal(predict(boot)$obj[[1]], rl(boot)[[1]],
              info=pst("predict.evmBoot: predict with type=rl gives same as direct call to rl with default arguments"))
  expect_equal(predict(boot, type="lp")$obj[[1]], linearPredictors(boot)[[1]],
              info=pst("predict.evmBoot: predict with type=lp gives same as direct call to linearPredictors with default arguments"))

  expect_equal(nm, length(pred$obj), info=pst("predict.evmBoot: output length"))
  expect_equal(paste("M.", from, sep=""), names(pred$obj)[1], info=pst("predict.evmBoot: names of output"))
  expect_equal(paste("M.", to, sep=""), names(pred$obj)[nm], info=pst("predict.evmBoot: names of output"))

  cnames <- c( "Mean","50%","2.5%","97.5%",names(X)[1:2])
  expect_equal(cnames, colnames(pred$obj[[1]]), info=pst("predict.evmBoot: colnames"))

  expect_equal(c(nx, 6), dim(pred$obj[[1]]), info=pst("predict.evmBoot: dimension"))
  for(i in 1:nm){
    expect_equal(unname(newX[, 1]), unname(pred$obj[[i]][,5]), info=pst("predict.evmBoot: covariates in output"))
    expect_equal(unname(newX[, 2]), unname(pred$obj[[i]][,6]), info=pst("predict.evmBoot: covariates in output"))
  }

  par(mfrow=n2mfrow(nx))
  plot(pred, sameAxes=FALSE, type="median", main="Bootstrap median rl")
  plot(pred, sameAxes=FALSE, type="mean", main="Bootstrap mean rl")
}
