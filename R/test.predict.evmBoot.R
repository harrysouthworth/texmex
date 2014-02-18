context("predict.evmBoot")

test_that("predict.evmBoot behaves as it should", {
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
    
  expect_that(target=predict(boot), equals(current=rl(boot)),   expect_that(target=predict(boot, equals(type="lp")),     
  expect_that(target=nm, equals(current=length(pred)),   expect_that(target=paste("M.", equals(from),   expect_that(target=paste("M.", equals(to),     
    cnames <- c( "Mean","50%","2.5%","97.5%",names(X)[1:2])
  expect_that(target=cnames, equals(current=colnames(pred[[1]])),     
  expect_that(target=c(nx, equals(6)),     for(i in 1:nm){
  expect_that(target=newX[, equals(1]),   expect_that(target=newX[, equals(2]),     }
    
    par(mfrow=n2mfrow(nx))
    plot(pred,sameAxes=FALSE,type="median",main="Bootstrap median rl")
    plot(pred,sameAxes=FALSE,type="mean",main="Bootstrap mean rl")
  }
}
)
