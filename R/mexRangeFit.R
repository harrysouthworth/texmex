`mexRangeFit` <-
function (x, which, quantiles=seq(0.5,0.9,length=9), start=c(.01, .01), R=10, nPass=3, trace=10,
          margins="laplace", constrain=TRUE, v=10)
{
  if (class(x) == "mex"){
    if( (!missing(margins))){
      warning("margins given, but already specified in 'mex' object.  Using 'mex' value")
    }
    if( (!missing(constrain))){
      warning("constrain given, but already specified in 'mex' object.  Using 'mex' value")
    }
    if( (!missing(v))){
      warning("v given, but already specified in 'mex' object.  Using 'mex' value")
    }
    if( (!missing(which))){
      warning("which given, but already specified in 'mex' object.  Using 'mex' value")
    }
    constrain <- x$dependence$constrain
    v <- x$dependence$v
    which <- x$dependence$which
    x <- x[[1]]
    margins <- x$margins
  } else {
    if (class(x) != "migpd"){
      stop("object should have class mex or migpd")
    }
    if (missing(which)) {
      which <- 1
      cat("Missing 'which'. Conditioning on", names(x$models)[which], ".\n")
    }
  }

  ests <- lapply(quantiles, function(qu, which, x, margins, start, constrain=constrain, v=v)
                                     mexDependence(x=x, which=which, dqu=qu, margins = margins, start=start, constrain=constrain, v=v),
                 which=which, x=x, margins = margins, start=start, constrain=constrain, v=v)

  boot <- lapply(ests, function(X, R, nPass, trace)
                                bootmex(x=X, R=R, nPass=nPass, trace=trace),
                 R=R, nPass=nPass, trace=trace)

  res <- list(ests=ests,boot=boot,quantiles=quantiles)
  oldClass(res) <- "mexRangeFit"
  res
}
 
print.mexRangeFit <- function(x, ...){
    out <- list(a = sapply(x$ests,function(x)x$dependence$coefficients[1,]),
                b = sapply(x$ests,function(x)x$dependence$coefficients[2,]))
    colnames(out$a) <- colnames(out$b) <- x$quantiles
    print(out)
    invisible()
}

summary.mexRangeFit <- function(object, ...){
  print(object)
}

plot.mexRangeFit <- function(x,col=2,bootcol="grey",addNexcesses=TRUE,...){
  ests <- x$ests
  boot <- x$boot
  quantiles <- x$quantiles
  PointEsts <- sapply(ests,function(X) coef(X$dependence))
  cof <- coef(ests[[1]]$dependence)
  whichName <- ests[[1]]$dependence$conditioningVariable
  which <- ests[[1]]$dependence$which
  data <- ests[[1]]$margins$data
  Names <- paste(rep(rownames(cof),dim(data)[2]-1),
                 paste(rep(colnames(cof),each=6),whichName,sep=" | "),sep="  ")
  R <- length(boot[[1]]$boot)
  
  for(i in 1:dim(PointEsts)[1]){
    if( sum((i %% 6) == 1:4) ){ # exclude plots from nuisance parameters m and s for which i mod 6 = 5,0 resp
      if(sum(PointEsts[i,])){
        Boot <- sapply(boot, function(x) sapply(x$boot, function(x) x$dependence[i]))
        ylim <- range(rbind(PointEsts[i,],Boot), na.rm=TRUE)
        plot(quantiles, PointEsts[i,], col=col, ylab=Names[i], type="b", ylim=ylim, ...)
        points(rep(quantiles,each=R),Boot,col=bootcol)
        points(quantiles, PointEsts[i,],col=col)
        if(addNexcesses){
          axis(3, at=axTicks(1), labels=sapply(axTicks(1), function(u,dat,which)sum(dat[,which] > quantile(dat[,which],u)),
                                               dat=data, which=which))
          mtext("# threshold excesses")
        }
      }
    }
  }
}

test.mexRangeFit <- function(){

  which <- 2
  quantiles <- seq(0.5,0.9,length=5)
  wmarmod <- migpd(winter, mqu=.7,  penalty="none")
  wmexmod.gum <- mexDependence(wmarmod, dqu=quantiles[1], margins="gumbel", constrain=FALSE,which=which)
  wmexmod.lap <- mexDependence(wmarmod, dqu=quantiles[1], margins="laplace",v=5,which=which)

  R <- 3
  mrf1 <- mexRangeFit(wmarmod,quantiles = quantiles,which=which,R=R,trace=R+1,v=5)
  mrf2 <- mexRangeFit(wmexmod.gum,quantiles = quantiles,R=R,trace=R+1)
  mrf3 <- mexRangeFit(wmexmod.lap,quantiles = quantiles,R=R,trace=R+1)

  checkException(mexRangeFit(TRUE,which=2),silent=TRUE,msg="mexRangeFit: exception handle")
  checkException(mexRangeFit(5,which=1),silent=TRUE,msg="mexRangeFit: exception handle")
  
  checkEquals(mrf1$ests[[1]][1:2],wmexmod.lap[1:2])
  checkEquals(mrf2$ests[[1]][1:2],wmexmod.gum[1:2])
  checkEquals(mrf3$ests[[1]][1:2],wmexmod.lap[1:2])
  
# now 2-d data

  mqu <- .7
  wavesurge.fit <- migpd(wavesurge,mqu=mqu)
  m <- mexDependence(wavesurge.fit,which=1,dqu=mqu)
  mrf4 <- mexRangeFit(wavesurge.fit,which=1,margins="laplace",R=R,trace=R+1)
  mrf5 <- mexRangeFit(m,R=R,trace=R+1)
  checkEquals(mrf4$ests[[2]][1:2],mrf5$ests[[2]][1:2])

# test specification of starting values
  R <- 5
  qu <- c(0.5,0.7,0.9)
  mrf6 <- mexRangeFit(wavesurge.fit,which=1,margins="laplace",constrain=TRUE, start=c(0.01,0.01),R=R,trace=R+1,quantiles = qu)
  mrf7 <- mexRangeFit(wavesurge.fit,which=2,margins="laplace",constrain=TRUE, start=m,R=R,trace=R+1,quantiles = qu)
  
# test plotting
  
  par(mfrow=c(2,2))
  plot(mrf6,main="start=(0.01,0.01)",addNexcesses=FALSE)
  plot(mrf7,main=paste("start=",signif(coef(m)$dependence[1:2],2)),addNexcesses=FALSE)
}
