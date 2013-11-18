test.predict.evmOpt <-
function(){
  
  for(Family in list(gpd,gev)){
    set.seed(20130513)
    pst <- function(msg) texmexPst(msg,Family=Family)
    
    # no covariates
    
    u    <- switch(Family$name,GPD=14,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    r.fit <- evm(data,th=u,family=Family)
    co <- coef(r.fit)
    
    if(Family$name == "GPD")checkEqualsNumeric(target=u,current = predict(r.fit,M=1/r.fit$rate)[[1]],msg=pst("predict.evmOpt: GPD retrieve threshold"))
    
    checkEquals(target=predict(r.fit), current=rl(r.fit),msg=pst("predict.evmOpt: predict with type=rl gives same as direct call to rl with default arguments"))
    checkEquals(target=predict(r.fit,type="lp"), current=linearPredictors(r.fit),msg=pst("predict.evmOpt: predict with type=lp gives same as direct call to linearpredictors with default arguments"))
    
    r.fit$rate <- 1
    prob <- c(0.5,0.9,0.95,0.99,0.999)
    for(p in prob){
      checkEqualsNumeric(target = Family$quant(p,t(co),r.fit),
                         current = unlist(predict(r.fit,M=1/(1-p))),msg=pst("predict.evmOpt: est ret levels no covariates"))
    }
    
    # with covariates
    
    n <- 1000
    M <- 1000
    
    # boundary cases with single covariates in only one parameter
    mu <- 1
    xi <- 0.05
    
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=cbind(X[,1],xi),GEV=cbind(mu,X[,1],xi))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,phi=~a,th=th,family=Family)
    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]
    pred <- predict(fit,M=M)
    
    checkEqualsNumeric(target=X$a,current=pred[[1]][,-1],msg=pst("predict.evmOpt: ret level correct reporting of covariates for single cov in phi only"))
    checkEqualsNumeric(target=Family$quant(1-1/M,AllCo,fit),
                       current = pred[[1]][,1],msg=pst("predict.evmOpt: ret level estimation with covariates in phi only"))
    
    mu <- 1
    sig <- 2
    
    param <- switch(Family$name,GPD=cbind(log(sig),X[,2]),GEV=cbind(mu,log(sig),X[,2]))
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,xi=~b,th=th,family=Family)
    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]
    pred <- predict(fit,M=M)
    
    checkEqualsNumeric(target=X$b,current=pred[[1]][,-1],msg=pst("predict.evmOpt: ret level correct reporting of covariates for single cov in xi only"))
    checkEqualsNumeric(target=Family$quant(1-1/M,AllCo,fit),
                       current = pred[[1]][,1],msg=pst("predict.evmOpt: ret level estimation with covariates in xi only"))
    
    # covariates in all parameters
    
    param <- switch(Family$name,GPD=cbind(X[,1],X[,2]),GEV=cbind(X[,1],X[,1],X[,2]))
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- switch(Family$name,
                  GPD = evm(Y, data=X,        phi=~a, xi=~b, th=th,family=Family),
                  GEV = evm(Y, data=X, mu=~a, phi=~a, xi=~b, th=th,family=Family))  
    
    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]
    
    checkEqualsNumeric(target=Family$quant(1-1/M,AllCo,fit),
                       current = predict(fit,M=M)[[1]][,1],msg=pst("predict.evmOpt: ret level estimation with covariates in all parameters"))
    
    # check multiple M
    M <- c(10,50,100,500,1000)
    
    target <- sapply(M,function(m,AllCo,fit) Family$quant(1-1/m,AllCo,fit),AllCo=AllCo,fit=fit)
    current <- predict(fit,M=M)
    
    for(i in 1:length(M)){
      checkEqualsNumeric(target[,i],current[[i]][,1],msg=pst("predict.evmOpt: ret level estimation multiple M"))
    }
    
    # new data
    nx <- 20
    M <- 1000
    newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
    AllCoNew <- predict(fit,type="lp",newdata=newX)$link[,1:length(Family$param)]
    
    checkEqualsNumeric(target=Family$quant(1-1/M,AllCoNew,fit),current=predict(fit,M=M,newdata=newX)[[1]][,1],msg=pst("predict.evmOpt: ret level ests with new data"))
    
    checkEqualsNumeric(dim(AllCoNew) + c(0,3),dim(predict(fit,ci=TRUE,newdata=newX)[[1]]), msg=pst("predict.evmOpt: dimension of return object for ci calc"))
    checkEqualsNumeric(dim(AllCoNew) + c(0,2),dim(predict(fit,se=TRUE,newdata=newX)[[1]]), msg=pst("predict.evmOpt: dimension of return object for se calc"))
    checkEqualsNumeric(dim(AllCoNew) + c(0,4),dim(predict(fit,se=TRUE,ci=TRUE,newdata=newX)[[1]]), msg=pst("predict.evmOpt: dimension of return object for se and ci calc"))
    
    Labels <- switch(Family$name,GPD=c("a","b"),GEV=c("a","a","b"))
    
    checkEquals(c("RL","2.5%","97.5%","se.fit",Labels), colnames(predict(fit,se=TRUE,ci=TRUE)[[1]]), msg=pst("predict.evmOpt: colnames of return obejct for se and ci calc, default alpha"))
    checkEquals(c("RL","5%","95%","se.fit",Labels), colnames(predict(fit,se=TRUE,ci=TRUE,alpha=0.1)[[1]]), msg=pst("predict.evmOpt: colnames of return obejct for se and ci calc, alpha=0.1"))
    
    # alpha
    
    alpha <- c(0.01,0.05,0.1,0.2,0.5,0.9,0.99)
    z <- matrix(qnorm(c(alpha/2,1-alpha/2)),ncol=2)
    
    for(a in 1:length(alpha)){
      p <- predict(fit,alpha=alpha[a],ci=TRUE,se=TRUE)[[1]]
      checkEquals(current = colnames(p)[2:3],target = c(paste(100*alpha[a]/2,"%",sep=""),paste(100*(1-alpha[a]/2),"%",sep="")),msg=pst("predict.evmOpt: labelling of confidence intervals"))
      checkEqualsNumeric(target = p[,1] + z[a,1]*p[,4],current = p[,2], msg=pst("predict.evmOpt: ret level Conf Interval calc for different alpha"))
      checkEqualsNumeric(target = p[,1] + z[a,2]*p[,4],current = p[,3], msg=pst("predict.evmOpt: ret level Conf Interval calc for different alpha"))
    }
    
    
    # linear predictors
    
    checkEqualsNumeric(target = c(nx,switch(Family$name,GPD=6,GEV=9)), dim(predict(fit,newdata=newX,se=TRUE,type="lp")$link), msg=pst("predict.evmOpt: dimension of return object, linear predictor, se calc"))
    checkEqualsNumeric(target = c(nx,switch(Family$name,GPD=8,GEV=12)),dim(predict(fit,newdata=newX,ci=TRUE,type="lp")$link), msg=pst("predict.evmOpt: dimension of return object, linear predictor, ci calc"))
    
    nameCI <- switch(Family$name,GPD = c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi","a","b"),GEV = c("mu","phi", "xi", "mu.lo","mu.hi","phi.lo", "phi.hi", "xi.lo", "xi.hi","a","a","b"))
    nameSE <- switch(Family$name,GPD = c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi","phi.se", "xi.se","a","b"),GEV = c("mu","phi", "xi", "mu.lo","mu.hi","phi.lo", "phi.hi", "xi.lo", "xi.hi","mu.se","phi.se", "xi.se","a","a","b"))
    checkEquals(target = nameCI, current = colnames(predict(fit,newdata=newX,ci=TRUE,type="lp")$link),msg=pst("predict.evmOpt: colnames for linear predictor return object"))
    checkEquals(target = nameSE, current = colnames(predict(fit,newdata=newX,ci=TRUE,se=TRUE,type="lp")$link),msg=pst("predict.evmOpt: colnames for linear predictor return object"))
    
    # unique
    
    newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1)) # checks for duplicates in one and both covariates.
    U <- !duplicated(newX)
    checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")$link,
                       target = predict(fit,newdata=newX,unique.=FALSE,type="lp")$link[U,], msg=pst("predict.evmOpt: functioning of argument unique, for linear predictor"))
    checkEqualsNumeric(current = predict(fit,newdata=newX)[[1]],
                       target =  predict(fit,newdata=newX,unique.=FALSE)[[1]][U,], msg=pst("predict.evmOpt: functioning of argument unique, for return levels"))
    
    # check standard errors - this takes a while since using bootstrap
    
    M <- c(10,100,500,1000,2000)
    newX <- data.frame("a"=rep(c(1,-1,2,-2),2),"b"=c(rep(0.1,4),rep(-0.1,4)))
    fit.p <- predict(fit, newdata=newX,se=TRUE,M=M)
    fit.seest <- unlist(lapply(fit.p,function(x) x[,2]))
    
    o <- options(warn=-1)
    fit.b <- evmBoot(fit,R=1000, trace=1100)
    options(o)
    fit.bp <- predict(fit.b,newdata=newX,all=TRUE,M=M)
    fit.seb <- lapply(fit.bp,function(X) apply(X,2,sd))
    fit.seboot <- unlist(fit.seb)
    
    checkTrue(all(abs((fit.seboot -  fit.seest) / fit.seest) < 0.3),msg=pst("predict.evmOpt: return level standard error estimate compared with bootstrap standard errors"))
  }  
}
