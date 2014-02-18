context("predict.evmSim")

test_that("predict.evmSim behaves as it should", {
    
  for(Family in list(gpd,gev)){
    
    pst <- function(msg) texmex:::texmexPst(msg,Family=Family)
    set.seed(20130513)
    # no covariates
    
    u    <- switch(Family$name,GPD=14,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    r.fit <- evm(data,th=u,family=Family,method="sim",trace=50000)
    co <- coef(r.fit)
    
    if(Family$name == "GPD"){
  expect_that(target=u, equals(current=predict(r.fit),     }
    
    r.fit$map$rate <- 1
    p <- c(0.5,0.9,0.95,0.99,0.999)
  expect_that(target=Family$quant(p, equals(t(co)),                        current = unlist(predict(r.fit,M=1/(1-p))),label=pst("predict.evmSim: ret level estimation"))
    
  expect_that(target=predict(r.fit, equals(M=1/(1-p))),   expect_that(target=predict(r.fit, equals(type="lp")$link),     
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
    
  expect_that(target=PostMeanRL(AllCo, equals(M)),     # multiple M
    
    M <- c(10,50,100,500,1000)
    
    target <- lapply(M,function(m) PostMeanRL(AllCo,m))
    current <- predict(fit,M=M)
    
    for(i in 1:length(M)){
  expect_that(target[[i]], equals(current[[i]][),     }
    
    # new data
    nx <- 20
    M <- 1000
    newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
    
    AllCoNew <- predict(fit,type="lp",newdata=newX,all=TRUE)$link
    
  expect_that(target=PostMeanRL(AllCoNew, equals(M)),   expect_that(target=as.matrix(newX), equals(current=predict(fit),     
    
    p <- predict(fit,all=TRUE,newdata=newX)
    alpha <- c(0.05,0.1)
    for(a in alpha){
      l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
      l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
      m.L <- unlist(l.L)
      m.U <- unlist(l.U)
      r <- predict(fit,newdata=newX,ci=TRUE,alpha=a)
  expect_that(target=m.L, equals(current=r[[1]][),   expect_that(target=m.U, equals(current=r[[1]][),     }
    
    # check linear predictors
    p <- predict(fit,type="lp",all=TRUE,newdata=newX)$link
    l <- lapply(p,function(l) apply(l,2,mean))
    m <- matrix(unlist(l),ncol=length(l[[1]]),byrow=TRUE)
    r <- predict(fit,type="lp",newdata=newX)$link
  expect_that(target=m, equals(current=r),     Offset <- switch(Family$name,GEV=1,GPD=0)
  expect_that(current=predict(fit, equals(newdata=newX),                        target = cbind(apply(cbind(rep(1,nx),newX[,1]) %*% t(fit$param[,Offset + 1:2]),1,mean),
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
  expect_that(target=m.L[, equals(1:npar]),   expect_that(target=m.U[, equals(1:npar]),     }
    
    # structure of output
    
  expect_that(target=c(n, equals(6)),     o <- options(warn=-1) # since se=TRUE gives a warning
  expect_that(target=c(n, equals(6)),     options(o)
    
  expect_that(target=c("Mean", equals("50%"),   expect_that(target=c("Mean", equals("50%"),     
  expect_that(target=c(nx, equals(4*npar+2)),     
    cnamesGPD <- c("phi: Mean", "phi: 50%", "phi: 2.5%", "phi: 97.5%", "xi: Mean", "xi: 50%", "xi: 2.5%", "xi: 97.5%")# this specific format assumed by plot.rl.evmSim and plot.lp.evmSim
    cnamesGEV <- c("mu: Mean",  "mu: 50%",  "mu: 2.5%",  "mu: 97.5%",  "phi: Mean", "phi: 50%", "phi: 2.5%", "phi: 97.5%", "xi: Mean", "xi: 50%", "xi: 2.5%", "xi: 97.5%")
    cnames<- switch(Family$name,GPD=cnamesGPD,GEV=cnamesGEV)
    
  expect_that(current=cnames, equals(target=colnames(predict(fit),     o <- options(warn=-1) # since se=TRUE gives a warning
  expect_that(current=cnames, equals(target=colnames(predict(fit),     options(o)
    
    # unique
    newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1))
    
  expect_that(current=predict(fit, equals(newdata=newX)[[1]]),   expect_that(current=predict(fit, equals(newdata=newX),   }
}
)
