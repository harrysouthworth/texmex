plot.lp.evmOpt <- function(x, main=NULL,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 1, polycol = 15, ...){
  family <- x$family
  x <- x$link
  
  if(dim(x)[1] == 1){
    stop("Need range of covariate values to plot linear predictors")
  }
  if(!any(colnames(x) == "phi.lo") ){
    stop("Please use ci.fit=TRUE in call to predict, to calculate confidence intervals")
  }

  Ests <- family$lp(data.frame(x))
  Names <- family$param
  cn <- colnames(x)
  which <- cn != "mu" & cn != "phi"    & cn != "xi" &
           cn != "mu.lo" & cn != "mu.hi" &
           cn != "phi.lo" & cn != "phi.hi" &
           cn != "xi.lo"  & cn != "xi.hi" &
           cn != "mu.se" & cn != "phi.se" & cn != "xi.se"

  X <- x[,which]
  if(is.null(dim(X))){
     X <- matrix(X)
     dimnames(X) <- list(dimnames(x)[[1]],dimnames(x)[[2]][which])
  }

  for(i in 1:length(Names)){
    for(j in 1:dim(X)[2]){
      if(length(unique(Ests[[i]][,1])) > 1){
        if(length(unique(X[,j])) > 1){
          ord <- order(X[,j])
          x <- X[ord,j]
          y <- Ests[[i]][ord,]
          plot(x, y[,1],type="n",ylab=Names[i],xlab=colnames(X)[j],main=main,ylim=range(y))

          if (polycol != 0){
            polygon(c( x,        rev(x)),
                    c(y[,2],rev(y[,3])),
                    col=polycol, border = FALSE) # Close polygon
          } else {
            lines(x, y[,2], col = cicol)
            lines(x, y[,3], col = cicol)
          }

          lines(x, y[,1], col = linecol[ 1 ] )
        }
      }
    }
  }
  invisible()
}

plot.lp.evmSim <- function(x, type="median", ...){
  if(dim(x$link)[1] == 1){
    stop("Need range of covariate values to plot linear predictors")
  }
  p <- x$family$param
  np <- length(p)
# re-format to same column structure as lp.evmOpt x
  ColIndexMeans <- 1+4*(0:(np-1))
  if(casefold(type) == "median"){
    offset <- 1
  } else if(casefold(type) == "mean") {
    offset <- 0
  } else {
    stop("type must be \"mean\" or \"median\" ")
  }
  which <- c(ColIndexMeans + offset,rep(1:2,np) + rep(ColIndexMeans+1,each=2), (4*np+1): dim(x$link)[2])
  x$link <- x$link[,which]
  colnames(x$link)[1:(3*np)] <-  c(p,paste(rep(p,each=2),rep(c(".lo",".hi"),np),sep=""))

  plot.lp.evmOpt(x,...)
}

plot.lp.evmBoot <- plot.lp.evmSim

test.plot.predict.evm <- function(){
# testing all of: plot.lp.evm* and plot.rl.evm* where * is opt, sim and boot

# first with no covariates
  n <- 100
  for(Family in list(gpd,gev)){
    set.seed(20130513)
    pst <- function(msg) texmexPst(msg,Family=Family)

    u    <- switch(Family$name,GPD=14,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    fit.opt <- evm(data,th=u,family=Family)
    fit.sim <- evm(data,th=u,method="sim",trace=100000,family=Family)
    fit.boot <- evmBoot(fit.opt,R=20,trace=30)
    
    M <- seq(20,1000,len=20)

    p.opt <- predict(fit.opt,M=M,ci=TRUE)
    p.sim <- predict(fit.sim,M=M,ci=TRUE)
    p.boot <- predict(fit.boot,M=M,ci=TRUE)

    par(mfrow=c(3,3))
    plot(p.opt,main=paste(Family$name,"\nMLE"))
    plot(p.sim,main=paste(Family$name,"\nMCMC"))
    plot(p.boot,main=paste(Family$name,"\nBootstrap"))

    p.lp.opt <- predict(fit.opt,type="lp",ci=TRUE)
    p.lp.sim <- predict(fit.sim,type="lp",ci=TRUE)
    p.lp.boot <- predict(fit.boot,type="lp",ci=TRUE)

    checkException(plot(p.lp.opt),silent=TRUE,msg=pst("plot.lp.evmOpt: fail if no covariates"))
    checkException(plot(p.lp.sim),silent=TRUE,msg=pst("plot.lp.evmSim: fail if no covariates"))
    checkException(plot(p.lp.boot),silent=TRUE,msg=pst("plot.lp.evmBoot: fail if no covariates"))

# now with covariates
    
    n <- 1000
    M <- 1000
    
    mu <- 1
    
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=cbind(X[,1],X[,2]),GEV=cbind(mu,X[,1],X[,2]))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    start <- switch(Family$name,GPD=c(0,1,0,1),GEV=c(mu,0,1,0,1))
    
    fit.opt <- evm(Y,data=X,phi=~a,xi=~b, th=th,family=Family,start=start)
    fit.sim <- evm(Y,data=X,phi=~a,xi=~b, th=th,family=Family,method="sim",trace=100000,start=start)
    o <- options(warn=-1)
    fit.boot <- evmBoot(fit.opt,R=20,trace=30)
    options(o)
    
    nx <- 3
    M <- seq(5,1000,len=20)
    newX <- data.frame(a=rnorm(nx),b=runif(nx,-0.1,0.1))

    p.opt <- predict(fit.opt,M=M,newdata=newX,ci=TRUE)
    p.sim <- predict(fit.sim,M=M,newdata=newX,ci=TRUE)
    p.boot <- predict(fit.boot,M=M,newdata=newX,ci=TRUE)

    p.lp.opt <- predict(fit.opt,type="lp",newdata=newX,ci=TRUE)
    p.lp.sim <- predict(fit.sim,type="lp",newdata=newX,ci=TRUE)
    p.lp.boot <- predict(fit.boot,type="lp",newdata=newX,ci=TRUE)

    par(mfrow=c(3,3))
    plot(p.opt,sameAxes=FALSE,main=paste(Family$name,"MLE\ndifferent axes"))
    plot(p.sim,sameAxes=FALSE,main=paste(Family$name,"MCMC\ndifferent axes"))
    plot(p.boot,sameAxes=FALSE,main=paste(Family$name,"Bootstrap\ndifferent axes"))

    plot(p.opt,sameAxes=TRUE,main=paste(Family$name,"MLE\nsame axes"))
    plot(p.sim,sameAxes=TRUE,main=paste(Family$name,"MCMC\nsame axes"))
    plot(p.boot,sameAxes=TRUE,main=paste(Family$name,"Bootstrap\nsame axes"))

    par(mfrow=c(4,4))
    plot(p.lp.opt,main=paste(Family$name,"MLE"))
    plot(p.lp.sim,main=paste(Family$name,"MCMC\nmean"),type="mean")
    plot(p.lp.sim,main=paste(Family$name,"MCMC\nmedian"),type="median")
    plot(p.lp.boot,main=paste(Family$name,"Bootstrap"))

# single covariate only:

    param <- switch(Family$name,GPD=cbind(X[1,1],X[,2]),GEV=cbind(mu,X[1,1],X[,2]))
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit.opt <- evm(Y,data=X,xi=~b,th=th,family=Family)
    fit.sim <- evm(Y,data=X,xi=~b,th=th,family=Family,method="sim",trace=100000)
    o <- options(warn=-1)
    fit.boot <- evmBoot(fit.opt,R=20,trace=30)
    options(o)

    p.opt <- predict(fit.opt,M=M,newdata=newX,ci=TRUE)
    p.sim <- predict(fit.sim,M=M,newdata=newX,ci=TRUE)
    p.boot <- predict(fit.boot,M=M,newdata=newX,ci=TRUE)

    p.lp.opt <- predict(fit.opt,type="lp",newdata=newX,ci=TRUE)
    p.lp.sim <- predict(fit.sim,type="lp",newdata=newX,ci=TRUE)
    p.lp.boot <- predict(fit.boot,type="lp",newdata=newX,ci=TRUE)

    par(mfrow=c(3,3))
    plot(p.opt,sameAxes=FALSE,main=paste(Family$name,"MLE"))
    plot(p.sim,sameAxes=FALSE,main=paste(Family$name,"MCMC"))
    plot(p.boot,sameAxes=FALSE,main=paste(Family$name,"Bootstrap"))

    par(mfrow=c(3,1))
    plot(p.lp.opt,main=paste(Family$name,"MLE"),polycol="cyan")
    plot(p.lp.sim,main=paste(Family$name,"MCMC"),polycol="cyan")
    plot(p.lp.boot,main=paste(Family$name,"Bootstrap"),polycol="cyan")
  }
}
