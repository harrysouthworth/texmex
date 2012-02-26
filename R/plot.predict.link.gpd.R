plot.lp.gpd <- function(x, main=NULL,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 1, polycol = 15, ...){

  if( !any(colnames(x) == "phi.lo") ){
    stop("Please use ci.fit=TRUE in call to predict, to calculate confidence intervals")
  }

  Ests <- list(phi=x[,c(1,3,4)],xi=x[,c(2,5,6)])
  Names <- c("phi","xi")
  cn <- colnames(x)
  which <- cn != "phi"    & cn != "xi" & 
           cn != "phi.lo" & cn != "phi.hi" & 
           cn != "xi.lo"  & cn != "xi.hi" & 
           cn != "phi.se" & cn != "xi.se"  
          
  X <- x[,which]
  if(is.null(dim(X))){
     X <- matrix(X)
     dimnames(X) <- list(dimnames(x)[[1]],dimnames(x)[[2]][which])
  }

  for(i in 1:2){
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

plot.lp.bgpd <- function(x, type="median", ...){
# re-format to same column structure as lp.gpd x
  if( casefold(type) == "median"){
    x <- x[,c(2,6,3,4,7,8,9:dim(x)[2])]
  } else if(casefold(type) == "mean") {
    x <- x[,c(1,5,3,4,7,8,9:dim(x)[2])]
  } else {
    stop("type must be \"mean\" or \"median\" ")
  }
  
  colnames(x)[1:6] <-  c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi")

  plot.lp.gpd(x,...)
}

plot.lp.bootgpd <- plot.lp.bgpd

test.plot.lp.gpd <- function(){
  n <- 100
  M <- 1000
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,exp(X[,1]),X[,2])
  X$Y <- Y
  fit <- gpd(Y,data=X,phi=~a, xi=~b,th=0)
  fitb <- gpd(Y,data=X,phi=~a, xi=~b,th=0,method="sim",trace=20000)
  o <- options(warn=-1)
  fit.boot <- bootgpd(fit,R=20,trace=30)
  options(o)
  
  nx <- 3
  M <- seq(5,1000,len=20)
  newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))

  p <- predict(fit,M=M,newdata=newX,ci=TRUE)
  pb <- predict(fitb,M=M,newdata=newX,ci=TRUE)
  pboot <- predict(fit.boot,M=M,newdata=newX,ci=TRUE)

  p.lp <- predict(fit,type="lp",newdata=newX,ci=TRUE)
  pb.lp <- predict(fitb,type="lp",newdata=newX,ci=TRUE)
  pboot.lp <- predict(fit.boot,type="lp",newdata=newX,ci=TRUE)
  
  par(mfrow=c(3,3))
  plot(p,sameAxes=FALSE,main="MLE")
  plot(pb,sameAxes=FALSE,main="MCMC")
  plot(pboot,sameAxes=FALSE,main="Bootstrap")
  
  par(mfrow=c(3,4))
  plot(p.lp,main="MLE")
  plot(pb.lp,main="MCMC")
  plot(pboot.lp,main="Bootstrap")
  
# single covariate only:

  Y <- rgpd(n,exp(X[1,1]),X[,2])
  X$Y <- Y
  fit <- gpd(Y,data=X,xi=~b,th=0)
  fitb <- gpd(Y,data=X,xi=~b,th=0,method="sim",trace=20000)
  o <- options(warn=-1)
  fit.boot <- bootgpd(fit,R=20,trace=30) 
  options(o)
  
  p <- predict(fit,M=M,newdata=newX,ci=TRUE)
  pb <- predict(fitb,M=M,newdata=newX,ci=TRUE)
  pboot <- predict(fit.boot,M=M,newdata=newX,ci=TRUE)
  
  p.lp <- predict(fit,type="lp",newdata=newX,ci=TRUE)
  pb.lp <- predict(fitb,type="lp",newdata=newX,ci=TRUE)
  pboot.lp <- predict(fit.boot,type="lp",newdata=newX,ci=TRUE)
  
  par(mfrow=c(3,3))
  plot(p,sameAxes=FALSE,main="MLE")
  plot(pb,sameAxes=FALSE,main="MCMC")
  plot(pboot,sameAxes=FALSE,main="Bootstrap")
  
  par(mfrow=c(3,1))
  plot(p.lp,main="MLE",polycol="cyan")
  plot(pb.lp,main="MCMC",polycol="cyan")
  plot(pboot.lp,main="Bootstrap",polycol="cyan")
}
