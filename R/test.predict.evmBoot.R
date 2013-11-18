test.predict.evmBoot <-
function(){
  # functionality all tested already in test.predict.evm, so just check output of correct format.
  
  n <- 1000
  nx <- 9
  R <- 10
  nm <- 20
  
  for(Family in list(gpd,gev)){
    
    pst <- function(msg) texmexPst(msg,Family=Family)
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
    
    checkEquals(target=predict(boot), current=rl(boot),msg=pst("predict.evmBoot: predict with type=rl gives same as direct call to rl with default arguments"))
    checkEquals(target=predict(boot,type="lp"), current=linearPredictors(boot),msg=pst("predict.evmBoot: predict with type=lp gives same as direct call to linearPredictors with default arguments"))
    
    checkEqualsNumeric(target=nm,current=length(pred),msg=pst("predict.evmBoot: output length"))
    checkEquals(target=paste("M.",from,sep=""),current=names(pred)[1],msg=pst("predict.evmBoot: names of output"))
    checkEquals(target=paste("M.",to,sep=""),current=names(pred)[nm],msg=pst("predict.evmBoot: names of output"))
    
    cnames <- c( "Mean","50%","2.5%","97.5%",names(X)[1:2])
    checkEquals(target=cnames,current=colnames(pred[[1]]),msg=pst("predict.evmBoot: colnames"))
    
    checkEqualsNumeric(target=c(nx,6),current=dim(pred[[1]]),msg=pst("predict.evmBoot: dimension"))
    for(i in 1:nm){
      checkEqualsNumeric(target=newX[,1],current=pred[[i]][,5],msg=pst("predict.evmBoot: covariates inoutput"))
      checkEqualsNumeric(target=newX[,2],current=pred[[i]][,6],msg=pst("predict.evmBoot: covariates inoutput"))
    }
    
    par(mfrow=n2mfrow(nx))
    plot(pred,sameAxes=FALSE,type="median",main="Bootstrap median rl")
    plot(pred,sameAxes=FALSE,type="mean",main="Bootstrap mean rl")
  }
}
