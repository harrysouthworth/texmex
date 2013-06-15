`mexRangeFit` <-
function (x, which, quantiles=seq(0.5,0.9,length=9), start=c(.01, .01), R=10, nPass=3, trace=10,
          col=2,bootcol="grey",margins="laplace", constrain=TRUE, v=10, addNexcesses=TRUE,...)
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

  PointEsts <- sapply(ests,function(X) coef(X$dependence))
  cof <- coef(ests[[1]]$dependence)
  whichName <- ests[[1]]$dependence$conditioningVariable
  Names <- paste(rep(rownames(cof),dim(x$data)[2]-1),
                 paste(rep(colnames(cof),each=6),whichName,sep=" | "),sep="  ")

  for(i in 1:dim(PointEsts)[1]){
    if( sum((i%%6) == 1:4) ){ # exclude plots from nuisance parameters m and s for which i mod 6 = 5,0 resp
      if(sum(PointEsts[i,])){
        Boot <- sapply(boot,function(x) sapply(x$boot,function(x) x$dependence[i]))
        ylim <- range(rbind(PointEsts[i,],Boot),na.rm=TRUE)
        plot(quantiles,PointEsts[i,],col=col,ylab=Names[i],type="b",ylim=ylim,...)
        points(rep(quantiles,each=R),Boot,col=bootcol)
        points(quantiles,PointEsts[i,],col=col)
        if(addNexcesses){
          axis(3,at=axTicks(1),labels=sapply(axTicks(1),function(u,dat,which)sum(dat[,which] > quantile(dat[,which],u)), dat=x$data,which=which))
          mtext("# threshold excesses")
        }
      }
    }
  }
}

test.mexRangeFit <- function(){

  wmarmod <- migpd(winter, mqu=.7,  penalty="none")
  wmexmod.gum <- mex(winter, mqu=.7,  penalty="none", margins="gumbel", constrain=FALSE, which=1)
  wmexmod.lap <- mex(winter, mqu=.7,  penalty="none", margins="laplace",v=5, which=1)

  par(mfrow=c(2,2))
  R <- 10
  mexRangeFit(wmarmod,which=1,margins="gumbel",constrain=FALSE,R=R,trace=R+1,
              main="Dependence threshold selection\nWinter data, Heffernan and Tawn 2004",cex=0.5,addNexcesses=FALSE)
  mexRangeFit(wmexmod.gum,main="Dependence threshold selection\nWinter data, Heffernan and Tawn 2004,\nGumbel margins",cex=0.5,addNexcesses=FALSE,R=R,trace=R+1)
  mexRangeFit(wmexmod.lap,main="Dependence threshold selection\nWinter data, Heffernan and Tawn 2004,\nLaplace margins",cex=0.5,addNexcesses=FALSE,R=R,trace=R+1)

  checkException(mexRangeFit(TRUE,which=2),silent=TRUE,msg="mexRangeFit: exception handle")
  checkException(mexRangeFit(5,which=1),silent=TRUE,msg="mexRangeFit: exception handle")

# now 2-d data

  wavesurge.fit <- migpd(wavesurge,mqu=.7)
  m <- mex(wavesurge,which=1,mqu=0.7)
  mexRangeFit(wavesurge.fit,which=1,margins="laplace",R=R,trace=R+1,
              main="Dependence threshold selection,\nwave and surge data, Coles 2001",addNexcesses=FALSE)
  mexRangeFit(wavesurge.fit,which=1,margins="gumbel",constrain=FALSE,R=R,trace=R+1,
              main="Dependence threshold selection,\nwave and surge data, Coles 2001",addNexcesses=FALSE)

# test specification of starting values
  R <- 5
  qu <- c(0.5,0.7,0.9)
  mexRangeFit(wavesurge.fit,which=1,margins="laplace",constrain=TRUE, start=c(0.01,0.01),R=R,trace=R+1,quantiles = qu,
              main="start=c(0.01,0.01)",addNexcesses=FALSE)
  mexRangeFit(wavesurge.fit,which=2,margins="laplace",constrain=TRUE, start=m,R=R,trace=R+1,quantiles = qu,
              main="start=fitted model",addNexcesses=FALSE)

}
