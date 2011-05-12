`mexRangeFit` <-
function (x, which, quantiles=seq(0.5,0.9,length=9), R=10, nPass=3, trace=10,
          col="red",bootcol="grey",margins="laplace", ...)
{

  if (class(x) == "mex"){
    x <- x[[2]]$migpd
    if( (!missing(margins))){
      warning("margins given, but already applied to 'mex' object.  Using 'mex' value")
    }
    margins <- x$margins
  } else if (class(x) != "migpd"){
    stop("object should have class mex or migpd")
  } 
  
  if (missing(which)) {
     cat("Missing 'which'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
     which <- 1
  }

  ests <- lapply(quantiles, function(qu, which, x, margins)
                                    mexDependence(x=x, which=which, dqu=qu, margins = margins),
                 which=which, x=x, margins = margins)
  boot <- lapply(quantiles, function(qu, x, which, R, nPass, trace, margins)
                                    bootmex(x=x, which=which, R=R, dqu=qu, nPass=nPass, trace=trace, margins=margins),
                 which=which, x=x, margins=margins, R=R, nPass=nPass, trace=trace)

  PointEsts <- sapply(ests,coef)
  cof <- coef(ests[[1]])
  whichName <- ests[[1]]$conditioningVariable  
  Names <- paste(rep(rownames(cof),dim(x$data)[2]-1),
                 paste(rep(colnames(cof),each=4),whichName,sep=" | "),sep="  ")

  for(i in 1:dim(PointEsts)[1]){
    if(sum(PointEsts[i,])){
      Boot <- sapply(boot,function(x) sapply(x$boot,function(x) x$dependence[i]))
      ylim <- range(rbind(PointEsts[i,],Boot),na.rm=TRUE)
      plot(quantiles,PointEsts[i,],col=col,ylab=Names[i],type="b",ylim=ylim,...)
      points(rep(quantiles,each=R),Boot,col=bootcol)
      points(quantiles,PointEsts[i,],col=col)
    }
  }
}

test.mexRangeFit <- function(){

  wmarmod <- migpd(winter, mqu=.7,  penalty="none")
  wmexmod.gum <- mex(winter, mqu=.7,  penalty="none", margins="gumbel")
  wmexmod.lap <- mex(winter, mqu=.7,  penalty="none", margins="laplace")
  
  par(mfrow=c(2,2))
  mexRangeFit(wmarmod,which=1,margins="gumbel",
              main="Dependence threshold selection\nWinter data, Heffernan and Tawn 2004",cex=0.5)
  mexRangeFit(wmexmod.gum,main="Dependence threshold selection\nWinter data, Heffernan and Tawn 2004,\nGumbel margins",cex=0.5)
  mexRangeFit(wmexmod.lap,main="Dependence threshold selection\nWinter data, Heffernan and Tawn 2004,\nLaplace margins",cex=0.5)
  
  o <- options()
  options(show.error.messages=FALSE)
  checkException(mexRangeFit(TRUE,which=2),msg="mexRangeFit: exception handle")
  checkException(mexRangeFit(5,which=1),msg="mexRangeFit: exception handle")
  options(op)
  
# now 2-d data

  wavesurge.fit <- migpd(wavesurge,mq=.7)
  mexRangeFit(wavesurge.fit,which=1,margins="laplace",
              main="Dependence threshold selection,\nwave and surge data, Coles 2001")
  mexRangeFit(wavesurge.fit,which=1,margins="gumbel",
              main="Dependence threshold selection,\nwave and surge data, Coles 2001")
  
  
}
