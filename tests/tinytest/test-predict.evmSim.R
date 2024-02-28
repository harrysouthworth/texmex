for(Family in list(gpd,gev)){

  pst <- function(msg) texmex:::texmexPst(msg,Family=Family)
  set.seed(20130513)
  # no covariates

  u    <- switch(Family$name,GPD=14,GEV=-Inf)
  data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)

  r.fit <- evm(data,th=u,family=Family,method="sim",trace=50000)
  co <- coef(r.fit)

  if(Family$name == "GPD"){
    expect_equal(u, predict(r.fit, M=1/r.fit$map$rate)$obj[[1]], info=pst("predict.evmSim: retrieve threshold"))
  }

  r.fit$map$rate <- 1
  p <- c(0.5,0.9,0.95,0.99,0.999)
  expect_equal(Family$quant(p, t(co), r.fit$map),
              unname(unlist(predict(r.fit,M=1/(1-p))$obj)), tolerance=0.01,
              info=pst("predict.evmSim: ret level estimation"))

  expect_equal(predict(r.fit, M=1/(1-p))$obj, unclass(rl(r.fit,M=1/(1-p))),
              info=pst("predict.evmSim: predict with type=rl gives same as direct call to rl with default arguments"))
  expect_equal(predict(r.fit, type="lp")$obj$link, linearPredictors(r.fit)$link,
              info=pst("predict.evmSim: predict with type=rl gives same as direct call to rl with default arguments"))

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

  AllCo <- predict(fit,type="lp",all=TRUE)$obj$link
  PostMeanRL <- function(AllCo,M){
    AllQuant <- lapply(AllCo,function(X)Family$quant(1-1/M,X[,1:length(Family$param)],fit$map))
    sapply(AllQuant,mean)
  }

  expect_equivalent(PostMeanRL(AllCo, M), predict(fit, M=M)$obj[[1]][,1],
              info=pst("predict.evmSim:retlevelestimationwithcovariatesinallparameters"))
  # multiple M

  M <- c(10,50,100,500,1000)

  target <- lapply(M,function(m) PostMeanRL(AllCo,m))
  current <- predict(fit,M=M)$obj

  for(i in 1:length(M)){
    expect_equal(target[[i]], unname(current[[i]][, 1]), tolerance=0.02,
                info=pst("predict.evmSim: ret level estimation multiple M"))
  }

  # new data
  nx <- 20
  M <- 1000
  newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))

  AllCoNew <- predict(fit,type="lp",newdata=newX,all=TRUE)$obj$link

  expect_equivalent(PostMeanRL(AllCoNew, M), predict(fit,M=M,newdata=newX)$obj[[1]][,1],
              info=pst("predict.evmSim: ret level ests with newdata"))
  expect_equal(c(as.matrix(newX)), c(predict(fit, M=M, newdata=newX)$obj[[1]][,2:3]),
              info=pst("predict.evmSim: ret level estimation with newdata, covariates added correctly to output"))


  p <- predict(fit,all=TRUE,newdata=newX)$obj
  alpha <- c(0.05,0.1)
  for(a in alpha){
    l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
    l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
    m.L <- unlist(l.L)
    m.U <- unlist(l.U)
    r <- predict(fit,newdata=newX,ci=TRUE,alpha=a)$obj
    expect_equal(unname(m.L), unname(r[[1]][, 3]), info=pst("predict.evmSim: lower conf ints for ret levels with newdata"))
    expect_equal(unname(m.U), unname(r[[1]][, 4]), info=pst("predict.evmSim: upper conf ints for ret levels with newdata"))
  }

  # check linear predictors
  p <- predict(fit,type="lp",all=TRUE,newdata=newX)$obj$link
  l <- lapply(p,function(l) apply(l,2,mean))
  m <- matrix(unlist(l),ncol=length(l[[1]]),byrow=TRUE)
  r <- predict(fit,type="lp",newdata=newX)$obj$link
  expect_equal(c(m), c(r), info=pst("predict.evmSim: linear predictors of parameters with newdata"))
  Offset <- switch(Family$name,GEV=1,GPD=0)
  expect_equal(unname(c(predict(fit, newdata=newX, type="lp")$obj$link[,Offset+1:2])),
              c(cbind(apply(cbind(rep(1,nx),newX[,1]) %*% t(fit$param[,Offset + 1:2]),1,mean),
                             apply(cbind(rep(1,nx),newX[,2]) %*% t(fit$param[,Offset + 3:4]),1,mean))),
              info = pst("predict.evmSim: linear predictor estimates"))

  alpha <- c(0.05,0.1)
  for(a in alpha){
    l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
    l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
    m.L <- matrix(unlist(l.L),ncol=length(l[[1]]),byrow=TRUE)
    m.U <- matrix(unlist(l.U),ncol=length(l[[1]]),byrow=TRUE)
    r <- predict(fit,type="lp",newdata=newX,ci=TRUE,alpha=a)$obj$link
    npar <- length(Family$param)
    expect_equal(c(m.L[, 1:npar]), c(r[,(1:(4*npar))[rep(c(F,F,T,F),npar)]]),
                info=pst("predict.evmSim: lower conf ints for linear predictors of parameters with newdata"))
    expect_equal(c(m.U[, 1:npar]), c(r[,(1:(4*npar))[rep(c(F,F,F,T),npar)]]),
                info=pst("predict.evmSim: upper conf ints for linear predictors of parameters with newdata"))
  }

  # structure of output

  expect_equal(c(n, 6), dim(predict(fit,ci=TRUE)$obj[[1]]),info=pst("predict.evmSim: dimension of output with ci calculation"))
  o <- options(warn=-1) # since se=TRUE gives a warning
  expect_equal(c(n, 6), dim(predict(fit,se=TRUE,ci=TRUE)$obj[[1]]),info=pst("predict.evmSim: dimension of output with ci and se calculation"))
  options(o)

  expect_equal(c("Mean", "50%", "2.5%","97.5%","a","b"), colnames(predict(fit,ci=TRUE)$obj[[1]]),
              info=pst("predict.evmSim: colnames of ret level ests with CI estimation"))
  expect_equal(c("Mean", "50%", "5%","95%","a","b"), colnames(predict(fit, alpha=0.1, ci=TRUE)$obj[[1]]),
              info=pst("predict.evmSim: colnames of ret level ests with CI estimation, alpha=0.1"))

  expect_equal(c(nx, 4*npar+2), dim(predict(fit,newdata=newX,ci=TRUE,type="lp")$obj$link),
              info=pst("predict.evmSim: dimension of linear predictor return object"))

  cnamesGPD <- c("phi:Mean", "phi:50%", "phi:2.5%", "phi:97.5%", "xi:Mean", "xi:50%", "xi:2.5%", "xi:97.5%")# this specific format assumed by plot.rl.evmSim and plot.lp.evmSim
  cnamesGEV <- c("mu:Mean",  "mu:50%",  "mu:2.5%",  "mu:97.5%",  "phi:Mean", "phi:50%", "phi:2.5%", "phi:97.5%", "xi:Mean", "xi:50%", "xi:2.5%", "xi:97.5%")
  cnames<- switch(Family$name,GPD=cnamesGPD,GEV=cnamesGEV)

  expect_equal(cnames, colnames(predict(fit, newdata=newX,ci=TRUE,type="lp")$obj$link)[1:(4*npar)],
              info=pst("predict.evmSim: colnames of lin predictors with CI calcs"))
  o <- options(warn=-1) # since se=TRUE gives a warning
  expect_equal(cnames, colnames(predict(fit, newdata=newX,ci=TRUE,se=TRUE,type="lp")$obj$link)[1:(4*npar)],
              info=pst("predict.evmSim: colnames of lin predictors with CI + SE calcs"))
  options(o)

  # unique
  newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1))

  expect_equal(predict(fit, newdata=newX)$obj[[1]], unique(predict(fit,newdata=newX,unique.=FALSE)$obj[[1]]),
              info=pst("predict.evmSim:uniquefunctioningforretlevelests"))
  expect_equal(predict(fit, newdata=newX, type="lp")$obj$link[,],
              unique(predict(fit,newdata=newX,unique.=FALSE,type="lp")$obj$link[,]),
              info=pst("predict.evmSim: unique functioning for lin pred ests"))
}
