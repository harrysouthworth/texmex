`mexRangeFit` <-
function (x, which, quantiles=seq(0.5,0.9,length=9), R=10, nPass=3, trace=10,
          col="red",bootcol="grey",...)
{

  if (class(x) == "mex"){
    x <- x[[1]]
  }

   if (missing(which)) {
       cat("Missing 'which'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
       which <- 1
   }

  ests <- lapply(quantiles,function(qu, which, x){
                                mexDependence(x=x, which=which, dqu=qu)
                           },
                           which=which, x=x)
  boot <- lapply(quantiles,function(qu, x, which, R, nPass, trace){
                                bootmex(x=x, which=which, R=R, dqu=qu, nPass=nPass, trace=trace)
                            },
                            x=x, which=which, R=R, nPass=nPass, trace=trace)

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
  
  par(mfrow=c(2,2))
  mexRangeFit(wmarmod,which=1,main="Dependence threshold selection\nWinter data, Heffernan and Tawn 2004",cex=0.5)
  
  checkException(mexRangeFit(TRUE,which=2),msg="mexRangeFit: exception handle")
  checkException(mexRangeFit(5,which=1),msg="mexRangeFit: exception handle")
  
# now 2-d data

  wavesurge.fit <- migpd(wavesurge,mq=.7)
  mexRangeFit(wavesurge.fit,which=1,main="Dependence threshold selection,\nwave and surge data, Coles 2001")
  
}
