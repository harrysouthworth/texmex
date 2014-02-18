context("predict.evmOpt")

test_that("predict.evmOpt behaves as it should", {
    
  for(Family in list(gpd,gev)){
    set.seed(20130513)
    pst <- function(msg) texmex:::texmexPst(msg,Family=Family)
    
    # no covariates
    
    u    <- switch(Family$name,GPD=14,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    r.fit <- evm(data,th=u,family=Family)
    co <- coef(r.fit)
    
    if(Family$name == "GPD")checkEqualsNumeric(target=u,current = predict(r.fit,M=1/r.fit$rate)[[1]],label=pst("predict.evmOpt: GPD retrieve threshold"))
    
  expect_that(target=predict(r.fit), equals(current=rl(r.fit)),   expect_that(target=predict(r.fit, equals(type="lp")),     
    r.fit$rate <- 1
    prob <- c(0.5,0.9,0.95,0.99,0.999)
    for(p in prob){
  expect_that(target=Family$quant(p, equals(t(co)),                          current = unlist(predict(r.fit,M=1/(1-p))),label=pst("predict.evmOpt: est ret levels no covariates"))
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
    
  expect_that(target=X$a, equals(current=pred[[1]][),   expect_that(target=Family$quant(1-1/M, equals(AllCo),                        current = pred[[1]][,1],label=pst("predict.evmOpt: ret level estimation with covariates in phi only"))
    
    mu <- 1
    sig <- 2
    
    param <- switch(Family$name,GPD=cbind(log(sig),X[,2]),GEV=cbind(mu,log(sig),X[,2]))
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,xi=~b,th=th,family=Family)
    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]
    pred <- predict(fit,M=M)
    
  expect_that(target=X$b, equals(current=pred[[1]][),   expect_that(target=Family$quant(1-1/M, equals(AllCo),                        current = pred[[1]][,1],label=pst("predict.evmOpt: ret level estimation with covariates in xi only"))
    
    # covariates in all parameters
    
    param <- switch(Family$name,GPD=cbind(X[,1],X[,2]),GEV=cbind(X[,1],X[,1],X[,2]))
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- switch(Family$name,
                  GPD = evm(Y, data=X,        phi=~a, xi=~b, th=th,family=Family),
                  GEV = evm(Y, data=X, mu=~a, phi=~a, xi=~b, th=th,family=Family))  
    
    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]
    
  expect_that(target=Family$quant(1-1/M, equals(AllCo),                        current = predict(fit,M=M)[[1]][,1],label=pst("predict.evmOpt: ret level estimation with covariates in all parameters"))
    
    # check multiple M
    M <- c(10,50,100,500,1000)
    
    target <- sapply(M,function(m,AllCo,fit) Family$quant(1-1/m,AllCo,fit),AllCo=AllCo,fit=fit)
    current <- predict(fit,M=M)
    
    for(i in 1:length(M)){
  expect_that(target[, equals(i]),     }
    
    # new data
    nx <- 20
    M <- 1000
    newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
    AllCoNew <- predict(fit,type="lp",newdata=newX)$link[,1:length(Family$param)]
    
  expect_that(target=Family$quant(1-1/M, equals(AllCoNew),     
  expect_that(dim(AllCoNew)+c(0, equals(3)),   expect_that(dim(AllCoNew)+c(0, equals(2)),   expect_that(dim(AllCoNew)+c(0, equals(4)),     
    Labels <- switch(Family$name,GPD=c("a","b"),GEV=c("a","a","b"))
    
  expect_that(c("RL", equals("2.5%"),   expect_that(c("RL", equals("5%"),     
    # alpha
    
    alpha <- c(0.01,0.05,0.1,0.2,0.5,0.9,0.99)
    z <- matrix(qnorm(c(alpha/2,1-alpha/2)),ncol=2)
    
    for(a in 1:length(alpha)){
      p <- predict(fit,alpha=alpha[a],ci=TRUE,se=TRUE)[[1]]
  expect_that(current=colnames(p)[2:3], equals(target=c(paste(100*alpha[a]/2),   expect_that(target=p[, equals(1]+z[a),   expect_that(target=p[, equals(1]+z[a),     }
    
    
    # linear predictors
    
  expect_that(target=c(nx, equals(switch(Family$name),   expect_that(target=c(nx, equals(switch(Family$name),     
    nameCI <- switch(Family$name,GPD = c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi","a","b"),GEV = c("mu","phi", "xi", "mu.lo","mu.hi","phi.lo", "phi.hi", "xi.lo", "xi.hi","a","a","b"))
    nameSE <- switch(Family$name,GPD = c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi","phi.se", "xi.se","a","b"),GEV = c("mu","phi", "xi", "mu.lo","mu.hi","phi.lo", "phi.hi", "xi.lo", "xi.hi","mu.se","phi.se", "xi.se","a","a","b"))
  expect_that(target=nameCI, equals(current=colnames(predict(fit),   expect_that(target=nameSE, equals(current=colnames(predict(fit),     
    # unique
    
    newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1)) # checks for duplicates in one and both covariates.
    U <- !duplicated(newX)
  expect_that(current=predict(fit, equals(newdata=newX),                        target = predict(fit,newdata=newX,unique.=FALSE,type="lp")$link[U,], label=pst("predict.evmOpt: functioning of argument unique, for linear predictor"))
  expect_that(current=predict(fit, equals(newdata=newX)[[1]]),                        target =  predict(fit,newdata=newX,unique.=FALSE)[[1]][U,], label=pst("predict.evmOpt: functioning of argument unique, for return levels"))
    
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
    
  expect_that(all(abs((fit.seboot-fit.seest)/fit.seest)<0.3), is_true(), abel=pst("predict.evmOpt:returnlevelstandarderrorestimatecomparedwithbootstrapstandarderrors"))
  }  
}
)
