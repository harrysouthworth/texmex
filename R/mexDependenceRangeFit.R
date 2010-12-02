`mexDependenceRangeFit` <-
function (x, which, quantiles=seq(0.5,0.9,length=9), nboot=10,col="red",bootcol="grey",...)
{
  ests <- lapply(quantiles,function(qu)mexDependence(x,which,gqu=qu))
  boot <- lapply(quantiles,function(qu)bootmex(x,which,R=nboot,gqu=qu))
  
  PointEsts <- sapply(ests,coef)
  cof <- coef(ests[[1]])
  whichName <- ests[[1]]$conditioningVariable  
  Names <- paste(rep(rownames(cof),dim(x$data)[2]-1),
                 paste(rep(colnames(cof),each=4),whichName,sep="|"),sep="  ")

  for(i in 1:dim(PointEsts)[1]){
    if(sum(PointEsts[i,])){
      Boot <- sapply(boot,function(x) sapply(x$boot,function(x) x$dependence[i]))
      ylim <- range(rbind(PointEsts[i,],Boot))
      plot(quantiles,PointEsts[i,],col=col,ylab=Names[i],type="b",ylim=ylim,...)
      points(rep(quantiles,each=nboot),Boot,col=bootcol)
      points(quantiles,PointEsts[i,],col=col)
    }
  }
}

test(mexDependenceRangeFit) <- function(){

  wmarmod <- migpd(winter, qu=.7,  penalty="none")
  
  par(mfrow=c(2,2))
  mexDependenceRangeFit(wmarmod,which=1,main="Dependence threshold selection\nWinter data, Heffernan and Tawn 2004",cex=0.5)
  
  checkException(mexDependenceRangeFit(TRUE,which=2),msg="mexDependenceRangeFit: exception handle")
  checkException(mexDependenceRangeFit(5,which=1),msg="mexDependenceRangeFit: exception handle")
}