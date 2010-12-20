'plot.mexDependence' <- function(x,quantiles=seq(0.1,by=0.2,len=5),col="grey",...){
   if ( class( x ) != "mexDependence" )
		 stop( "you need to use an object with class 'mexDependence'" )
   
   z <- x$Z
   n <- dim(z)[1]

   for(i in 1:(dim(z)[2])){
      p <- seq(x$dqu,1-1/n,length=n)
      plot(p,z[,i],xlab=paste("F(",x$conditioningVariable,")",sep=""),
           ylab=paste("Z   ",colnames(z)[i]," | ",x$conditioningVariable,sep=""),col=col,...)
      lines(lowess(p,z[,i]),col="red")
      plot(p,abs(z[,i] - mean(z[,i])),xlab=paste("F(",x$conditioningVariable,")",sep=""),
           ylab=paste("|Z - mean(Z)|   ",colnames(z)[i]," | ",x$conditioningVariable,sep=""),col=col,...)
      lines(lowess(p,abs(z[,i] - mean(z[,i]))),col="blue")
   }

   xmax <- max(x$migpd$data[,x$which])
   sig <- coef(x$migpd)[3,x$which]
   xi <- coef(x$migpd)[4,x$which]
   marThr <- x$migpd$mth[x$which]
   marP   <- x$migpd$mqu[x$which]
   if(xi < 0) upperEnd <- marThr - sig/xi
   len <- 1001

  for(i in 1:(dim(z)[2])){
      depThr <- c(quantile(x$migpd$data[,x$which],x$dqu))
      d <- xmax-depThr
      xlim <- c(depThr-0.1*d, depThr + 1.5*d)
      names(xlim) <- NULL

      SetPlim <- TRUE
      if(xi < 0 && xlim[2] > upperEnd){
        xlim[2] <-  upperEnd
        plim <- 1
        SetPlim <- FALSE
      }
      
      if( SetPlim ) plim <- pgpd(xlim[2],sigma=sig,xi=xi,u=marThr)
      
      plotp <- seq(x$dqu,plim,len=len)[-len] # take out largest point to avoid Inf in next line if xlim[2]==upperEnd
      co <- coef(x)[,i]
      xq <- -log(-log(plotp))
      zq <- quantile(x$Z[,i],quantiles)
      yq <- sapply(zq, function(z)co["a"] * xq + co["c"] - co["d"]*log(xq) + xq^co["b"] * z)
      plotx <- revGumbel(plotp,data=x$migpd$data[,x$which],qu=marP,th=marThr,sigma=sig,xi=xi)
      ploty <- apply(exp(-exp(-yq)),2,revGumbel,data=as.matrix(x$migpd$data[,-x$which])[,i],
                     qu=x$migpd$mqu[-x$which][i],th=x$migpd$mth[-x$which][i],
                     sigma=coef(x$migpd)[3,-x$which][i],xi=coef(x$migpd)[4,-x$which][i])
      
      plot(x$migpd$data[,x$which],as.matrix(x$migpd$data[,-x$which])[,i], xlab=x$conditioningVariable,ylab=colnames(z)[i],col=col,...)
      abline(v=depThr)
      for(j in 1:length(quantiles)){
        lines(plotx,ploty[,j],lty=2)
      }
    }
}

  
test.plot.mexDependence <- function(){

smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
wmarmod <- migpd(winter, mqu=.7,  penalty="none")

mySdep1 <- mexDependence(smarmod,which=1, dqu=0.7)
myWdep1 <- mexDependence(wmarmod,which=1, dqu=0.75)
mySdep2 <- mexDependence(smarmod,which=2, dqu=0.8)
myWdep2 <- mexDependence(wmarmod,which=2, dqu=0.85)
mySdep3 <- mexDependence(smarmod,which=3, dqu=0.9)
myWdep3 <- mexDependence(wmarmod,which=3, dqu=0.7)
mySdep4 <- mexDependence(smarmod,which=4, dqu=0.75)
myWdep4 <- mexDependence(wmarmod,which=4, dqu=0.8)
mySdep5 <- mexDependence(smarmod,which=5, dqu=0.85)
myWdep5 <- mexDependence(wmarmod,which=5, dqu=0.9)

# check plots produced for variety of thresholds and parameter combinations
par(mfrow=c(3,4),pty="m")
plot(mySdep1,main=paste("Summer\nfitting quantile =",mySdep1$dqu))
plot(myWdep1,main=paste("Winter\nfitting quantile =",myWdep1$dqu))
plot(mySdep2,main=paste("Summer\nfitting quantile =",mySdep2$dqu))
plot(myWdep2,main=paste("Winter\nfitting quantile =",myWdep2$dqu))
plot(mySdep3,main=paste("Summer\nfitting quantile =",mySdep3$dqu))
plot(myWdep3,main=paste("Winter\nfitting quantile =",myWdep3$dqu))
plot(mySdep4,main=paste("Summer\nfitting quantile =",mySdep4$dqu))
plot(myWdep4,main=paste("Winter\nfitting quantile =",myWdep4$dqu))
plot(mySdep5,main=paste("Summer\nfitting quantile =",mySdep5$dqu))
plot(myWdep5,main=paste("Winter\nfitting quantile =",myWdep5$dqu))
   
   checkException(plot.mexDependence(smarmod),msg="mexDependence: exception handle")
   checkException(plot.mexDependence(TRUE),msg="mexDependence: exception handle")
   
# check execution for 2-d data
wavesurge.fit <- migpd(wavesurge,mqu=0.8)
wavesurge.mex <- mexDependence(wavesurge.fit,dqu=0.8,which=2)
par(mfrow=c(2,2))
res <- plot(wavesurge.mex,main="Wave surge data")
checkEquals(res,NULL,msg = "plot.mexDependence: successful execution")
}
