extremalIndex <- function(y,data=NULL,threshold)
# intevals estimator of the Extremal Index, Ferro and Segers JRSS B (2003)
# assumes data points equally spaced in time and no missing data (ie missing time points)
{
  if (!missing(data)) {
     y <- ifelse(deparse(substitute(y))== "substitute(y)", deparse(y),deparse(substitute(y)))
     y <- formula(paste(y, "~ 1"))
     y <- model.response(model.frame(y, data=data))
  }
  theCall <- match.call()
  timeAll <- 1:length(y)
  thExceedance <- y > threshold
  thExceedanceProb <- mean(thExceedance)
  nExceed <- sum(thExceedance)
  exceedanceTimes <- timeAll[thExceedance]
  interExceedTimes <- exceedanceTimes[-1] - exceedanceTimes[-nExceed]

  if ( max(interExceedTimes) > 2 ) {
    Est <- (2 * sum(interExceedTimes - 1)^2) / ( (nExceed-1) * sum( (interExceedTimes - 1) * (interExceedTimes-2) ) )
  } else {
    Est <- (2 * sum(interExceedTimes)^2) / ( (nExceed-1) * sum( interExceedTimes^2) )
  }

  EIintervals <- min(1,Est)

  res <- list(EIintervals = EIintervals,
              threshold = threshold,
              TotalN = length(y),
              nExceed = nExceed,
              thExceedanceProb = thExceedanceProb,
              call = theCall,
              interExceedTimes = interExceedTimes,
              thExceedances = y[thExceedance],
              exceedanceTimes = timeAll[thExceedance],
              y = y,
              data = data)

  oldClass(res) <- "extremalIndex"

  res
}

print.extremalIndex <- function(x,...)
{
  cat("\nLength of original series",x$TotalN,"\n")
  cat("Threshold", x$threshold,"\n")
  cat("Number of Threshold Exceedances",x$nExceed,"\n")
  cat("Intervals estimator of Extremal Index", x$EIintervals,"\n")
}

show.extremalIndex <- print.extremalIndex

plot.extremalIndex <- function(x,...)
{
  NormInterExceedTimes <- x$interExceedTimes * x$thExceedanceProb

  StdExpQuantiles <- qexp(ppoints(NormInterExceedTimes))
  Theta <- x$EIintervals

  plot(StdExpQuantiles, sort(NormInterExceedTimes),xlab="Standard Exponential Quantiles",ylab="Interexceedance Times",cex=0.7,...)
  abline(v=qexp(1-Theta))
  abline(a = -qexp(1-Theta)/Theta, b=1/Theta)
  title(paste("Threshold=",x$threshold))
}

declust <- function(y, r=NULL, data=NULL, ...)
{
  if (!missing(data)) {
     y <- deparse(substitute(y))
     y <- formula(paste(y, "~ 1"))
     y <- model.response(model.frame(y, data=data))
  }
  UseMethod("declust",y)
}

declust.default <- function(y,r=NULL,data=NULL,verbose=TRUE,...)
{
  if(missing(data)){
    ei <- extremalIndex(y,...)
  } else {
    ei <- extremalIndex(substitute(y),data,...)
  }
  if(verbose & is.null(r)){
    print(ei)
  }

  declust(ei,r=r)
}

declust.extremalIndex <- function(y,r=NULL,...)
{
  theCall <- match.call()
  Times <- y$interExceedTimes
  sTimes <- sort(Times, decreasing=TRUE)

  if(is.null(r)){
    C <- floor(y$EIintervals * y$nExceed) + 1
    C <- min(C,length(Times)) # needed if short series and C < number of interexceedance times
    while(sTimes[C-1] == sTimes[C]) C <- C-1
    r <- sTimes[C]
    method <- "intervals"
  } else {
    method <- "runs"
  }

  clusters <- rep(1,length(y$thExceedances))
  clusters[-1] <- 1+cumsum(Times > r)
  sizes <- tabulate(clusters)
  C <- max(clusters)

  clusterMaxima <- sapply(1:C,function(i) max(y$thExceedances[clusters == i]))
  isClusterMax <- rep(FALSE,length(clusters))
  for(i in 1:C){
    isClusterMax[clusters == i & y$thExceedances == max(y$thExceedances[clusters == i])][1] <- TRUE
  }

  res <- list(clusters = clusters,
              sizes=sizes,
              clusterMaxima = clusterMaxima,
              isClusterMax = isClusterMax,
              y = y$y,
              data = y$data,
              threshold=y$threshold,
              EIintervals = y$EIintervals,
              call=theCall,
              InterExceedTimes=Times,
              InterCluster = Times > sTimes[C],
              thExceedances = y$thExceedances,
              exceedanceTimes = y$exceedanceTimes,
              r=r, nClusters = C, method=method)

  oldClass(res) <- "declustered"

  res
}


print.declustered <- function(x,...){
  print(x$call)
  cat("\nThreshold ",x$threshold,"\n")
  cat("Declustering using the",x$method,"method, run length",x$r,"\n")
  cat("Identified",length(x$sizes),"clusters.\n")
}

show.declustered <- print.declustered

plot.declustered <- function(x,ylab="Data",...){
  plot(x$y,xlab="",ylab=ylab)
  abline(h=x$threshold,col=2)
  for(i in 1:length(x$sizes)){
    points(x$exceedanceTimes[x$clusters == i],x$thExceedances[x$clusters == i],col=2,type="b")
  }
}

bootExtremalIndex <- function(x){
  if( class(x) == "extremalIndex"){
    x <- declust(x)
  } else if(class(x) != "declustered"){
    stop("x must be of class extremalIndex or declust")
  }
  nc <- length(x$sizes)
  boot.interExceedTimes <- boot.thExceedances <- NULL
  for(i in 1:nc){
    clust <- sample(unique(x$clusters),1)
    boot.interExceedTimes <- c(boot.interExceedTimes,diff(x$exceedanceTimes[x$clusters == clust]))
    boot.thExceedances <- c(boot.thExceedances,x$thExceedances[x$clusters == clust])
    if(i < nc){
      boot.interExceedTimes <- c(boot.interExceedTimes,sample(x$InterExceedTimes[x$InterCluster], 1))
    }
  }

  boot.exceedanceTimes <- cumsum(c(1,boot.interExceedTimes))
  boot.data <- rep(-1,max(boot.exceedanceTimes))
  boot.data[boot.exceedanceTimes] <- boot.thExceedances

  boot.data
}

extremalIndexRangeFit <- function(y,data=NULL,umin=quantile(y,.5),umax=quantile(y,0.95),nint=10,nboot=100,alpha=.05,xlab="Threshold",addNexcesses=TRUE, estGPD=TRUE, verbose=TRUE, trace=10, ...){

  if (!missing(data)) {
     y <- deparse(substitute(y))
     y <- formula(paste(y, "~ 1"))
     y <- model.response(model.frame(y, data=data))
  }

  EI <- SH <- SC <-
    list(m=numeric(nint),boot=matrix(0,nrow=nint,ncol=nboot))

  u <- seq(umin, umax, length = nint)
  for (i in 1:nint) {
    if(verbose){
      cat("\n", i,"of",nint,"thresholds: bootstrap ... ")
    }
    z <- extremalIndex(y,threshold=u[i])
    EI$m[i] <- z$EIintervals
    d <- declust(z)
    if(estGPD){
      gpd.d <- gpd.declustered(d)
      co.d <- coef(gpd.d)
      SH$m[i] <- co.d[2]
      SC$m[i] <- exp(co.d[1]) - co.d[2]*u[i]
    }

    for(j in 1:nboot){
      if(verbose & j %% trace == 0){
        cat(j,"")
      }
      boot <- bootExtremalIndex(d)
      z.b <- extremalIndex(boot,threshold=u[i])
      EI$boot[i,j] <- z.b$EIintervals
      if(estGPD){
        z.d <- declust(z.b)
        z.d$clusterMaxima <- rgpd(z.d$nClusters,exp(co.d[1]),co.d[2],u=z.d$threshold)
        gpd.b <- try(gpd.declustered(z.d,cov="numeric"))
        if(class(gpd.b) == "try-error"){
          SH$boot[i,j] <- SC$boot[i,j] <- NA
        } else {
          SH$boot[i,j] <- coef(gpd.b)[2]
          SC$boot[i,j] <- exp(coef(gpd.b)[1]) -  SH$boot[i,j]*u[i]
        }
      }
    }
  }
  EI$ul <- apply(EI$boot,1,quantile,alpha/2)
  EI$up <- apply(EI$boot,1,quantile,1-alpha/2)
  if(estGPD){
    SC$ul <- apply(SC$boot,1,quantile,alpha/2,na.rm=TRUE)
    SC$up <- apply(SC$boot,1,quantile,1-alpha/2,na.rm=TRUE)
    SH$ul <- apply(SH$boot,1,quantile,alpha/2,na.rm=TRUE)
    SH$up <- apply(SH$boot,1,quantile,1-alpha/2,na.rm=TRUE)
  }

  plots <- function(l,...){
    plot(u, l$m, ylim=c(min(l$ul),max(l$up)),type = "b", xlab=xlab, ...)
    for (i in 1:nint) lines(c(u[i], u[i]), c(l$ul[i], l$up[i]))
    if(addNexcesses){
      axis(3,at=axTicks(1),labels=sapply(axTicks(1),function(u)max(declust(extremalIndex(y,threshold=u))$clusters)))
      mtext("# threshold excesses")
    }
  }

  plots(EI,main="Extremal Index",ylab=expression(theta),...)
  if(estGPD){
    plots(SC,main="Scale parameter",ylab=expression(sigma),...)
    plots(SH,main="Shape parameter",ylab=expression(xi),...)
  }

  invisible()
}

gpd.declustered <- function(y, ...){
  theCall <- match.call()
  if(is.null(y$data)){
    res <- gpd(y$clusterMaxima, th = y$threshold, ...)
  } else {
    response <- y$clusterMaxima
    dat <- cbind(response,y$data[y$y>y$threshold,][y$isClusterMax,])
    res <- gpd(response, data=dat, th = y$threshold, ...)
  }

  clusterRate <- max(y$clusters) / length(y$y)
  if(class(res) == "gpd"){
    res$rate <- clusterRate
  } else if(class(res) == "bgpd") {
    res$map$rate <- clusterRate
  }
  res$call <- theCall
  res
}

test.extremalIndex <- function(){
  tol <- 0.0001
  th <- quantile(rain,seq(0.7,0.99,len=10))
  for(i in 1:length(th)){
    texmex.ei <- extremalIndex(rain,threshold=th[i])
    Ferro.ei  <- texmex:::.extRemes.exi.intervals(rain > th[i])

    Ferro.clust <- texmex:::.extRemes.decluster.intervals(rain> th[i], Ferro.ei)
    texmex.clust <- declust(texmex.ei)

    Ferro.runs <-  texmex:::.extRemes.decluster.runs(rain> th[i], 3)
    texmex.runs <- declust(rain,threshold=th[i],r=3,verbose=FALSE)

    checkEqualsNumeric(texmex.ei$EIintervals, Ferro.ei,
          tolerance = tol,msg="extremalIndex: extRemes implementation")
    checkEqualsNumeric(texmex.clust$sizes, Ferro.clust$size,
          tolerance = tol,msg="extremalIndex: declustering")

    checkEqualsNumeric(texmex.runs$nCluster,Ferro.runs$nc,msg="extremalIndex: runs declustering nc")
    checkEqualsNumeric(texmex.runs$sizes,Ferro.runs$size,msg="extremalIndex: runs declustering sizes")
    }

# check passing data through data frames

  data <- data.frame(RAIN=rain[1:1000], rnorm=rnorm(1000), num=1:1000)
  extremalIndexRangeFit(RAIN, data,verbose=FALSE,nboot=20,nint=7)
  extremalIndexRangeFit(data$RAIN,verbose=FALSE,nboot=20,nint=7)

  data.de <- declust(RAIN,data=data,th=th[1],verb=FALSE)
  resp.de <- declust(data$RAIN,th=th[1],verb=FALSE)

  data.ei <- extremalIndex(RAIN,data=data,threshold=th[1])
  resp.ei <- extremalIndex(data$RAIN,threshold=th[1])

  checkEqualsNumeric(data.ei$EIintervals,resp.ei$EIintervals,tolerance=tol,msg="extremalIndex: using data frame to pass response")
  checkEqualsNumeric(data.de$clusters,resp.de$clusters,tolerance=tol,msg="extremalIndex: using data frame to pass numeric response to declustering")

}
