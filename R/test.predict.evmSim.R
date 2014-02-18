test.predict.evmSim <-
function(){
  
  for(Family in list(gpd,gev)){
    
    pst <- function(msg) texmexPst(msg,Family=Family)
    set.seed(20130513)
    # no covariates
    
    u    <- switch(Family$name,GPD=14,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    r.fit <- evm(data,th=u,family=Family,method="sim",trace=50000)
    co <- coef(r.fit)
    
    if(Family$name == "GPD"){
      checkEqualsNumeric(target=u,current=predict(r.fit,M=1/r.fit$map$rate)[[1]], msg=pst("predict.evmSim: retrieve threshold"))
    }
    
    r.fit$map$rate <- 1
    p <- c(0.5,0.9,0.95,0.99,0.999)
    checkEqualsNumeric(target = Family$quant(p,t(co),r.fit$map), tolerance=0.01,
                       current = unlist(predict(r.fit,M=1/(1-p))),msg=pst("predict.evmSim: ret level estimation"))
    
    checkEquals(target=predict(r.fit,M=1/(1-p)), current=rl(r.fit,M=1/(1-p)),msg=pst("predict.evmSim: predict with type=rl gives same as direct call to rl with default arguments"))
    checkEquals(target=predict(r.fit,type="lp")$link, current=linearPredictors(r.fit)$link,msg=pst("predict.evmSim: predict with type=rl gives same as direct call to rl with default arguments"))
    
    # with covariates
    
    n <- 1000
    M <- 1000
    mu <- 1
    
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=X,GEV=cbind(mu,X))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    start <- switch(Family$name,GPD=c(0,1,0,1),GEV=c(1,0,1,0,1))
    fit <- evm(Y,data=X,phi=~a,xi=~b,th=th,method="sim",trace=50000,family=Family,start=start)
    
    AllCo <- predict(fit,type="lp",all=TRUE)$link
    PostMeanRL <- function(AllCo,M){
      AllQuant <- lapply(AllCo,function(X)Family$quant(1-1/M,X[,1:length(Family$param)],fit$map))
      sapply(AllQuant,mean)
    }
    
    checkEqualsNumeric(target=PostMeanRL(AllCo,M),current = predict(fit,M=M)[[1]][,1],msg=pst("predict.evmSim: ret level estimation with covariates in all parameters"))
    # multiple M
    
    M <- c(10,50,100,500,1000)
    
    target <- lapply(M,function(m) PostMeanRL(AllCo,m))
    current <- predict(fit,M=M)
    
    for(i in 1:length(M)){
      checkEqualsNumeric(target[[i]],current[[i]][,1],tolerance=0.02,msg=pst("predict.evmSim: ret level estimation multiple M"))
    }
    
    # new data
    nx <- 20
    M <- 1000
    newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
    
    AllCoNew <- predict(fit,type="lp",newdata=newX,all=TRUE)$link
    
    checkEqualsNumeric(target=PostMeanRL(AllCoNew,M),current=predict(fit,M=M,newdata=newX)[[1]][,1],msg=pst("predict.evmSim: ret level ests with new data"))
    checkEqualsNumeric(target = as.matrix(newX), current = predict(fit,M=M,newdata=newX)[[1]][,2:3],msg=pst("predict.evmSim: ret level estimation with new data, covariates added correctly to output"))
    
    
    p <- predict(fit,all=TRUE,newdata=newX)
    alpha <- c(0.05,0.1)
    for(a in alpha){
      l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
      l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
      m.L <- unlist(l.L)
      m.U <- unlist(l.U)
      r <- predict(fit,newdata=newX,ci=TRUE,alpha=a)
      checkEqualsNumeric(target = m.L,current = r[[1]][,3],msg=pst("predict.evmSim : lower conf ints for ret levels with new data"))
      checkEqualsNumeric(target = m.U,current = r[[1]][,4],msg=pst("predict.evmSim : upper conf ints for ret levels with new data"))
    }
    
    # check linear predictors
    p <- predict(fit,type="lp",all=TRUE,newdata=newX)$link
    l <- lapply(p,function(l) apply(l,2,mean))
    m <- matrix(unlist(l),ncol=length(l[[1]]),byrow=TRUE)
    r <- predict(fit,type="lp",newdata=newX)$link
    checkEqualsNumeric(target = m,current = r,msg=pst("predict.evmSim : linear predictors of parameters with new data"))
    Offset <- switch(Family$name,GEV=1,GPD=0)
    checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")$link[,Offset +1:2],
                       target = cbind(apply(cbind(rep(1,nx),newX[,1]) %*% t(fit$param[,Offset + 1:2]),1,mean),
                                      apply(cbind(rep(1,nx),newX[,2]) %*% t(fit$param[,Offset + 3:4]),1,mean)),
                       msg = pst("predict.evmSim: linear predictor estimates"))
    
    alpha <- c(0.05,0.1)
    for(a in alpha){
      l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
      l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
      m.L <- matrix(unlist(l.L),ncol=length(l[[1]]),byrow=TRUE)
      m.U <- matrix(unlist(l.U),ncol=length(l[[1]]),byrow=TRUE)
      r <- predict(fit,type="lp",newdata=newX,ci=TRUE,alpha=a)$link
      npar <- length(Family$param)
      checkEqualsNumeric(target = m.L[,1:npar],current = r[,(1:(4*npar))[rep(c(F,F,T,F),npar)]],msg=pst("predict.evmSim : lower conf ints for linear predictors of parameters with new data"))
      checkEqualsNumeric(target = m.U[,1:npar],current = r[,(1:(4*npar))[rep(c(F,F,F,T),npar)]],msg=pst("predict.evmSim : upper conf ints for linear predictors of parameters with new data"))
    }
    
    # structure of output
    
    checkEqualsNumeric(target = c(n,6), current = dim(predict(fit,ci=TRUE)[[1]]), msg=pst("predict.evmSim: dimension of output with ci calculation"))
    o <- options(warn=-1) # since se=TRUE gives a warning
    checkEqualsNumeric(target = c(n,6), current = dim(predict(fit,se=TRUE,ci=TRUE)[[1]]), msg=pst("predict.evmSim: dimension of output with ci and se calculation"))
    options(o)
    
    checkEquals(target = c("Mean", "50%","2.5%","97.5%","a","b"), colnames(predict(fit,ci=TRUE)[[1]]), msg=pst("predict.evmSim: colnames of ret level ests with CI estimation"))
    checkEquals(target = c("Mean", "50%","5%","95%","a","b"), colnames(predict(fit,alpha=0.1,ci=TRUE)[[1]]), msg=pst("predict.evmSim: colnames of ret level ests with CI estimation, alpha=0.1"))
    
    checkEqualsNumeric(target = c(nx,4*npar+2), dim(predict(fit,newdata=newX,ci=TRUE,type="lp")$link), msg=pst("predict.evmSim: dimension of linear predictor return object"))
    
    cnamesGPD <- c("phi: Mean", "phi: 50%", "phi: 2.5%", "phi: 97.5%", "xi: Mean", "xi: 50%", "xi: 2.5%", "xi: 97.5%")# this specific format assumed by plot.rl.evmSim and plot.lp.evmSim
    cnamesGEV <- c("mu: Mean",  "mu: 50%",  "mu: 2.5%",  "mu: 97.5%",  "phi: Mean", "phi: 50%", "phi: 2.5%", "phi: 97.5%", "xi: Mean", "xi: 50%", "xi: 2.5%", "xi: 97.5%")
    cnames<- switch(Family$name,GPD=cnamesGPD,GEV=cnamesGEV)
    
    checkEquals(current = cnames, target = colnames(predict(fit,newdata=newX,ci=TRUE,type="lp")$link)[1:(4*npar)], msg=pst("predict.evmSim: col names of lin predictors with CI calcs"))
    o <- options(warn=-1) # since se=TRUE gives a warning
    checkEquals(current = cnames, target = colnames(predict(fit,newdata=newX,ci=TRUE,se=TRUE,type="lp")$link)[1:(4*npar)], msg=pst("predict.evmSim: col names of lin predictors with CI+SE calcs"))
    options(o)
    
    # unique
    newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1))
    
    checkEqualsNumeric(current = predict(fit,newdata=newX)[[1]], target = unique(predict(fit,newdata=newX,unique.=FALSE)[[1]]),pst("predict.evmSim: unique functioning for ret level ests"))
    checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")$link[,], target = unique(predict(fit,newdata=newX,unique.=FALSE,type="lp")$link[,]),msg=pst("predict.evmSim: unique functioning for lin pred ests"))
  }
}
