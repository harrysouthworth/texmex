plotrl.gpd <- # intended as a diagnostic for a gpd fitted with no covariates. Called by plot.gpd
function(object, alpha = .050,
         xlab, ylab, main,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 0, polycol = 15, smooth = TRUE,RetPeriodRange=NULL ){

    if(dim(object$X.phi)[2] > 1 | dim(object$X.xi)[2] > 1){
      stop("use plot method for object returned by predict.gpd to see rl plots if covariates in model")
    }
    if (missing(xlab) || is.null(xlab)) { xlab <- "Return period" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Return level" }
    if (missing(main) || is.null(main)) { main <- "Return Level Plot" }

    xdat <- object$y
    n <- length(xdat) / object$rate # Number of obs prior to thresholding

    if(!is.null(RetPeriodRange)){
      ran <- log10(RetPeriodRange)
      jj <- seq(ran[1],ran[2],by=0.1)
    } else {
      jj <- seq(-1, max(3.75,log10(n)),by=0.1)
    }

    m <- unique( c(1/object$rate, 10^jj) )
    xm <- matrix(unlist(rl(object,M=m,ci.fit=TRUE,alpha=alpha)),ncol=3,byrow=TRUE)

    U <- object$threshold - abs(object$threshold/100)
    plotX <- xm[,1] > U

    xrange <- range(m)
    yrange <- range(c(xdat, range(xm[plotX,])))

    plotRLgpd(m[plotX],xm[plotX,],polycol,cicol,linecol,ptcol,n,xdat,pch,
              smooth,xlab,ylab,main,xrange=xrange,yrange=yrange)

    invisible(list(m=m, xm=xm))
}

plot.rl.gpd <- function(x, # method for rl.(boot or b)gpd object, which may have covariates.  Plots return level for each unique row in design matrix
         xlab, ylab, main,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 0, polycol = 15, smooth = TRUE, sameAxes=TRUE, type="median", ...){
    if (missing(xlab) || is.null(xlab)) { xlab <- "Return period" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Return level" }
    if (missing(main) || is.null(main)) {
      main <- "Return Level Plot"
      SetMain <- TRUE
    } else {
      if(length(main) == 1){
        SetMain <- TRUE
      } else {
        SetMain <- FALSE
      }
    }

    nm <- length(names(x))
    if(nm < 2) {
      stop("Need to have more than one value of M at which to plot return level curve")
    }
    nd <- dim(x[[1]])[2]
    ncov <- length(unlist(x)) / (nm * nd)
    ValNames <- colnames(x[[1]])

    if(!SetMain & length(main) != ncov){
      stop("main must be length 1 or number of unique covariates for prediction")
    }

    if(length(ValNames) == 1 |  length(grep("%",ValNames)) < 2){
      stop("Please use ci.fit=TRUE in call to predict, to calculate confidence intervals")
    }

    Array <- array(unlist(x),c(ncov,nd,nm),dimnames=list(NULL,ValNames,names(x)))

    if( class(x) == "rl.gpd"){
      if(any(dimnames(x[[1]])[[2]] == "se.fit")){
        which <- dimnames(x[[1]])[[2]] != "se.fit"
        nd <- nd-1
        ValNames <- ValNames[which]
        Unlist <- unlist(x)[rep(which,each=ncov)]
      } else {
        Unlist <- unlist(x)
      }
    } else if(class(x) == "rl.bgpd" | class(x) == "rl.bootgpd"){
        if(casefold(type) == "median"){
          which <- dimnames(x[[1]])[[2]] != "Mean"
        } else if(casefold(type) == "mean") {
          which <- dimnames(x[[1]])[[2]] != "50%"
        } else {
          stop("type must be \"mean\" or \"median\" ")
        }
        nd <- nd-1
        ValNames <- ValNames[which]
        Unlist <- unlist(x)[rep(which,each=ncov)]
    }

    Array <- array(Unlist,c(ncov,nd,nm),dimnames=list(NULL,ValNames,names(x)))

    m <- as.numeric(substring(dimnames(Array)[[3]],first=3))

    if(ncov>1){
      covnames <- dimnames(Array)[[2]][-(1:3)]
    } else {
      covnames <- ""
    }

    if(sameAxes){
       yrange <- range(Array[,1:3,])
    }

    for(i in 1:ncov){
      xm <- t(Array[i,,])
      cov <- xm[1,-(1:3)]

      if(!sameAxes){
        yrange <- range(xm[,1:3])
      }

      if(SetMain){
        if(length(covnames) == 1 && covnames != ""){
          Main <- paste(main,"\n", paste(covnames,"=",signif(cov,2),collapse=", "))
        } else {
          Main <- main
        }
      } else {
        Main <- main[i]
      }
      plotRLgpd(m,xm,polycol = polycol,cicol=cicol,linecol=linecol,ptcol=ptcol,pch=pch,
                smooth=smooth,xlab=xlab,ylab=ylab,main=Main,xrange=range(m),yrange=yrange)
    }

    invisible(list(m=m,xm=Array))
}

plot.rl.bootgpd <- plot.rl.bgpd <- plot.rl.gpd

plotRLgpd <- function(M,xm,polycol,cicol,linecol,ptcol,n,xdat,pch,smooth,xlab,ylab,main,xrange,yrange){
# worker function - called by plotrl.gpd, plot.rl.gpd, plot.rl.bgpd

    o <- order(M) # in case the return period are not in ascending order.
    M <- M[o]
    xm <- xm[o,]

    plot(M, xm[,1], log = "x", type = "n",
         xlim=xrange, ylim=yrange, xlab = xlab, ylab = ylab, main = main)

      if (smooth) {
        splo <- spline(log(M), xm[,2] , 200)
        sphi <- spline(log(M), xm[,3] , 200)
        if ( polycol != 0 ) {
            polygon( exp(c( splo$x, rev( sphi$x ) )),
	                       c( splo$y, rev( sphi$y ) ),
                     col = polycol , border=FALSE)
         } # Close if (polycol
         lines( exp(splo$x), splo$y, col = cicol )
         lines( exp(sphi$x), sphi$y, col = cicol )
      } else{
        if (polycol != 0){
            polygon(c( M,        rev( M)),
                    c(xm[,2],rev(xm[,3])),
                    col=polycol, border = FALSE) # Close polygon
        } else {
            lines(M, xm[,2], col = cicol)
            lines(M, xm[,3], col = cicol)
        }
      }

      lines(M, xm[,1], col = linecol[ 1 ] )

    # Add observed data to the plot
    if(!missing(xdat) & !missing(n)){
      ly <- length(xdat)
      points(1 / (1 - ((n - ly + 1):n) / (n + 1)), sort(xdat), pch=pch, col=ptcol)
      box()
    }
 }

test.plotrl.gpd <- function()
{
# no covariates

  rain.fit <- gpd(rain,th=14)
  par(mfrow=c(1,1))
  plotrl.gpd(rain.fit)

# with covariates

  n <- 100
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,exp(X[,1]),X[,2])
  X$Y <- Y
  fit <- gpd(Y,data=X,phi=~a,xi=~b,th=0)
  rl <- predict(fit,ci=TRUE)

  checkException(plotrl.gpd(fit),silent=TRUE,msg="plotrl.gpd : failure for model with covariates")
  checkException(plot(rl),silent=TRUE,msg="plot.rl.gpd : failure when use only 1 value of M for RL calc")

  nx <- 6
  M <- seq(5,500,length=30)
  newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
  rl <- predict(fit,newdata=newX,ci.fit=TRUE,se=TRUE,M=M)
  par(mfrow=n2mfrow(nx))
  plot(rl,sameAxes=TRUE,main=paste("Validation suite plot",1:nx),polycol="cyan")
  plot(rl,sameAxes=FALSE,polycol="magenta")

  checkException(plot(predict(fit,newdata=newX,ci=FALSE)),silent=TRUE,msg="plotrl.gpd: failure if no conf ints supplied")

}

