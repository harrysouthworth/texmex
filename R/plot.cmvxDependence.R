'plot.mexDependence' <- function(x,quantiles=seq(0.1,by=0.2,len=5),col="grey",...){
   if ( class( x ) != "mexDependence" )
		 stop( "you need to use an object with class 'mexDependence'" )
   
   z <- x$Z
   n <- dim(z)[1]
   for(i in 1:(dim(z)[2])){
      p <- seq(x$gqu[i],1-1/n,length=n)
      plot(p,z[,i],xlab=paste("F(",x$conditioningVariable,")",sep=""),
           ylab=paste("Z   ",colnames(z)[i],"|",x$conditioningVariable,sep=""),col=col,...)
    abline(h=c(-2,0,2),lty=2)
   }

   xmax <- max(x$migpd$data[,x$which])
   sig <- coef(x$migpd)[3,x$which]
   xi <- coef(x$migpd)[4,x$which]
   marThr <- x$migpd$th[x$which]
   marP   <- x$migpd$qu[x$which]
   if(xi < 0) upperEnd <- marThr - sig/xi
   len <- 1001

  for(i in 1:(dim(z)[2])){
      depThr <- c(quantile(x$migpd$data[,x$which],x$gqu[i]))
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

      plotp <- seq(x$gqu[i],plim,len=len)[-len] # take out largest point to avoid Inf in next line if xlim[2]==upperEnd
      co <- coef(x)[,i]
      xq <- -log(-log(plotp))
      zq <- quantile(x$Z[,i],quantiles)
      yq <- sapply(zq, function(z)co["a"] * xq + co["c"] - co["d"]*log(xq) + xq^co["b"] * z)
      plotx <- revGumbel(plotp,data=x$migpd$data[,x$which],qu=marP,th=marThr,sigma=sig,xi=xi)
      ploty <- apply(exp(-exp(-yq)),2,revGumbel,data=x$migpd$data[,-x$which][,i],
                     qu=x$migpd$qu[-x$which][i],th=x$migpd$th[-x$which][i],
                     sigma=coef(x$migpd)[3,-x$which][i],xi=coef(x$migpd)[4,-x$which][i])
      
      plot(x$migpd$data[,x$which],x$migpd$data[,-x$which][,i], xlab=x$conditioningVariable,ylab=colnames(z)[i],col=col,...)
      abline(v=depThr)
      for(j in 1:length(quantiles)){
        lines(plotx,ploty[,j],lty=2)
      }
    }
}

test(plot.mexDependence) <- function(){

smarmod <- migpd(summer, qu=c(.9, .7, .7, .85, .7), penalty="none")
wmarmod <- migpd(winter, qu=.7,  penalty="none")

mySdep1 <- mexDependence(smarmod,which=1, gqu=0.7)
myWdep1 <- mexDependence(wmarmod,which=1, gqu=0.8)
mySdep2 <- mexDependence(smarmod,which=2, gqu=0.9)
myWdep2 <- mexDependence(wmarmod,which=2, gqu=0.75)
mySdep3 <- mexDependence(smarmod,which=3, gqu=0.75)
myWdep3 <- mexDependence(wmarmod,which=3, gqu=0.85)
mySdep4 <- mexDependence(smarmod,which=4, gqu=0.85)
myWdep4 <- mexDependence(wmarmod,which=4, gqu=0.73)
mySdep5 <- mexDependence(smarmod,which=5, gqu=0.73)
myWdep5 <- mexDependence(wmarmod,which=5, gqu=0.71)

# check plots produced for variety of thresholds and parameter combinations
par(mfrow=c(2,2))
plot(mySdep1,main="Summer")
plot(myWdep1,main="Winter")
plot(mySdep2,main="Summer")
plot(myWdep2,main="Winter")
plot(mySdep3,main="Summer")
plot(myWdep3,main="Winter")
plot(mySdep4,main="Summer")
plot(myWdep4,main="Winter")
plot(mySdep5,main="Summer")
plot(myWdep5,main="Winter")
   
   checkException(plot.mexDependence(smarmod))
   checkException(plot.mexDependence(TRUE))
 
}