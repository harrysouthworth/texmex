'plot.mex' <- function(x,quantiles=seq(0.1,by=0.2,len=5),col="grey",...){
   if ( class( x ) != "mex" )
		 stop( "you need to use an object with class 'mex'" )
   
   mar <- x[[1]]
   dep <- x[[2]]
   z <- dep$Z
   n <- dim(z)[1]
   xmax <- max(mar$data[,dep$which])
   sig <- coef(mar)[3,dep$which]
   xi <- coef(mar)[4,dep$which]
   marThr <- mar$mth[dep$which]
   marP   <- mar$mqu[dep$which]
   if(xi < 0) upperEnd <- marThr - sig/xi
   len <- 1001
   
   for(i in 1:(dim(z)[2])){
      p <- seq(dep$dqu,1-1/n,length=n)
      plot(p,z[,i],xlab=paste("F(",dep$conditioningVariable,")",sep=""),
           ylab=paste("Z   ",colnames(z)[i]," | ",dep$conditioningVariable,sep=""),col=col,...)
      lines(lowess(p,z[,i]),col="red")
      plot(p,abs(z[,i] - mean(z[,i])),xlab=paste("F(",dep$conditioningVariable,")",sep=""),
           ylab=paste("|Z - mean(Z)|   ",colnames(z)[i]," | ",dep$conditioningVariable,sep=""),col=col,...)
      lines(lowess(p,abs(z[,i] - mean(z[,i]))),col="blue")

      depThr <- c(quantile(mar$data[,dep$which],dep$dqu))
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
      
      plotp <- seq(dep$dqu,plim,len=len)[-len] # take out largest point to avoid Inf in next line if xlim[2]==upperEnd
      co <- coef(dep)[,i]
      xq <- -log(-log(plotp))
      zq <- quantile(dep$Z[,i],quantiles)
      yq <- sapply(zq, function(z, co, xq){
      						co["a"] * xq + co["c"] - co["d"]*log(xq) + xq^co["b"] * z
      					}, # Close function
      					 xq=xq, co=co
      			   ) # Close sapply
      
      plotx <- revTransform(plotp,data=mar$data[,dep$which],qu=marP,th=marThr,sigma=sig,xi=xi)
      ploty <- apply(exp(-exp(-yq)),2,revTransform,data=as.matrix(mar$data[,-dep$which])[,i],
                     qu=mar$mqu[-dep$which][i],th=mar$mth[-dep$which][i],
                     sigma=coef(mar)[3,-dep$which][i],xi=coef(mar)[4,-dep$which][i])
      
      plot(mar$data[,dep$which],as.matrix(mar$data[,-dep$which])[,i], xlab=dep$conditioningVariable,ylab=colnames(z)[i],col=col,...)
      abline(v=depThr)
      for(j in 1:length(quantiles)){
        lines(plotx,ploty[,j],lty=2)
      }
    }
}

  
test.plot.mex <- function(){

smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
wmarmod <- migpd(winter, mqu=.7,  penalty="none")

mySdep1.lap <- mexDependence(smarmod,which=1, dqu=0.7)
myWdep1.lap <- mexDependence(wmarmod,which=1, dqu=0.75)
mySdep2.lap <- mexDependence(smarmod,which=2, dqu=0.8)
myWdep2.lap <- mexDependence(wmarmod,which=2, dqu=0.85)
mySdep3.lap <- mexDependence(smarmod,which=3, dqu=0.9)
myWdep3.lap <- mexDependence(wmarmod,which=3, dqu=0.7)
mySdep4.lap <- mexDependence(smarmod,which=4, dqu=0.75)
myWdep4.lap <- mexDependence(wmarmod,which=4, dqu=0.8)
mySdep5.lap <- mexDependence(smarmod,which=5, dqu=0.85)
myWdep5.lap <- mexDependence(wmarmod,which=5, dqu=0.9)

mySdep1.gum <- mexDependence(smarmod,which=1, dqu=0.7, margins = "gumbel",constrain=FALSE)
myWdep1.gum <- mexDependence(wmarmod,which=1, dqu=0.75, margins = "gumbel",constrain=FALSE)
mySdep2.gum <- mexDependence(smarmod,which=2, dqu=0.8, margins = "gumbel",constrain=FALSE)
myWdep2.gum <- mexDependence(wmarmod,which=2, dqu=0.85, margins = "gumbel",constrain=FALSE)
mySdep3.gum <- mexDependence(smarmod,which=3, dqu=0.9, margins = "gumbel",constrain=FALSE)
myWdep3.gum <- mexDependence(wmarmod,which=3, dqu=0.7, margins = "gumbel",constrain=FALSE)
mySdep4.gum <- mexDependence(smarmod,which=4, dqu=0.75, margins = "gumbel",constrain=FALSE)
myWdep4.gum <- mexDependence(wmarmod,which=4, dqu=0.8, margins = "gumbel",constrain=FALSE)
mySdep5.gum <- mexDependence(smarmod,which=5, dqu=0.85, margins = "gumbel",constrain=FALSE)
myWdep5.gum <- mexDependence(wmarmod,which=5, dqu=0.9, margins = "gumbel",constrain=FALSE)


# check plots produced for variety of thresholds and parameter combinations
par(mfcol=c(3,4),pty="m")
plot(mySdep1.gum,main=paste("Summer\nfitting quantile =",mySdep1.gum$dependence$dqu,"\nGumbel margins"))
plot(myWdep1.gum,main=paste("Winter\nfitting quantile =",myWdep1.gum$dependence$dqu,"\nGumbel margins"))
plot(mySdep2.gum,main=paste("Summer\nfitting quantile =",mySdep2.gum$dependence$dqu,"\nGumbel margins"))
plot(myWdep2.gum,main=paste("Winter\nfitting quantile =",myWdep2.gum$dependence$dqu,"\nGumbel margins"))
plot(mySdep3.gum,main=paste("Summer\nfitting quantile =",mySdep3.gum$dependence$dqu,"\nGumbel margins"))
plot(myWdep3.gum,main=paste("Winter\nfitting quantile =",myWdep3.gum$dependence$dqu,"\nGumbel margins"))
plot(mySdep4.gum,main=paste("Summer\nfitting quantile =",mySdep4.gum$dependence$dqu,"\nGumbel margins"))
plot(myWdep4.gum,main=paste("Winter\nfitting quantile =",myWdep4.gum$dependence$dqu,"\nGumbel margins"))
plot(mySdep5.gum,main=paste("Summer\nfitting quantile =",mySdep5.gum$dependence$dqu,"\nGumbel margins"))
plot(myWdep5.gum,main=paste("Winter\nfitting quantile =",myWdep5.gum$dependence$dqu,"\nGumbel margins"))
   
plot(mySdep1.lap,main=paste("Summer\nfitting quantile =",mySdep1.lap$dependence$dqu,"\nLaplace margins"))
plot(myWdep1.lap,main=paste("Winter\nfitting quantile =",myWdep1.lap$dependence$dqu,"\nLaplace margins"))
plot(mySdep2.lap,main=paste("Summer\nfitting quantile =",mySdep2.lap$dependence$dqu,"\nLaplace margins"))
plot(myWdep2.lap,main=paste("Winter\nfitting quantile =",myWdep2.lap$dependence$dqu,"\nLaplace margins"))
plot(mySdep3.lap,main=paste("Summer\nfitting quantile =",mySdep3.lap$dependence$dqu,"\nLaplace margins"))
plot(myWdep3.lap,main=paste("Winter\nfitting quantile =",myWdep3.lap$dependence$dqu,"\nLaplace margins"))
plot(mySdep4.lap,main=paste("Summer\nfitting quantile =",mySdep4.lap$dependence$dqu,"\nLaplace margins"))
plot(myWdep4.lap,main=paste("Winter\nfitting quantile =",myWdep4.lap$dependence$dqu,"\nLaplace margins"))
plot(mySdep5.lap,main=paste("Summer\nfitting quantile =",mySdep5.lap$dependence$dqu,"\nLaplace margins"))
plot(myWdep5.lap,main=paste("Winter\nfitting quantile =",myWdep5.lap$dependence$dqu,"\nLaplace margins"))
   
op <- options()
options(show.error.messages=FALSE)
   checkException(plot.mex(smarmod),msg="plot.mex: exception handle")
   checkException(plot.mex(TRUE),msg="plot.mex: exception handle")
options(op)
   
# check execution for 2-d data
wavesurge.fit <- migpd(wavesurge,mqu=0.8)
wavesurge.mex <- mexDependence(wavesurge.fit,dqu=0.8,which=2,margins="gumbel",constrain=FALSE)
par(mfrow=c(2,2))
res <- plot(wavesurge.mex,main="Wave surge data")
checkEquals(res,NULL,msg = "plot.mex: successful execution")
}
