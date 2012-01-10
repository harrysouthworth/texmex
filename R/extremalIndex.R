extremalIndex <- function(data,threshold)
# intevals estimator of the Extremal Index, Ferro and Segers JRSS B (2003)
# assumes data points equally spaced in time and no missing data (ie missing time points)
{
  theCall <- match.call()
  timeAll <- 1:length(data)
  thExceedance <- data > threshold
  thExceedanceProb <- mean(thExceedance)
  nExceed <- sum(thExceedance)
  exceedanceTimes <- timeAll[thExceedance]
  interExceedTimes <- exceedanceTimes[-1] - exceedanceTimes[-nExceed]
  
  if ( max(interExceedTimes) > 2 ) {
    Est <- (2 * sum(interExceedTimes - 1)^2) / ( (nExceed-1) * sum( (interExceedTimes - 1) * (interExceedTimes-2) ) )
  } else {
    Est <- (2 * sum(interExceedTimes)^2) / ( (nExceed-1) * sum( interExceedTimes^2) )
  }
  
  EstIntervals <- min(1,Est)
  
  res <- list(EstIntervals = EstIntervals,
              threshold = threshold,
              TotalN = length(data),
              nExceed = nExceed,
              thExceedanceProb = thExceedanceProb,
              call = theCall,
              interExceedTimes = interExceedTimes,
              thExceedances = data[thExceedance],
              exceedanceTimes = timeAll[thExceedance],
              data = data)
  
  oldClass(res) <- "extremalIndex"
  
  res
}

print.extremalIndex <- function(Est)
{
  print(Est$call) 
  cat("\nLength of original series",Est$TotalN,"\n")
  cat("Threshold", Est$threshold,"\n")
  cat("Number of Threshold Exceedances",Est$nExceed,"\n")
  cat("Intervals estimator of Extremal Index", Est$EstIntervals,"\n")
}

show.extremalIndex <- print.extremalIndex 

plot.extremalIndex <- function(Est)
{
  NormInterExceedTimes <- Est$interExceedTimes * Est$thExceedanceProb
  
  StdExpQuantiles <- qexp(ppoints(NormInterExceedTimes))
  Theta <- Est$EstIntervals
  
  plot(StdExpQuantiles, sort(NormInterExceedTimes),xlab="Standard Exponential Quauntiles",ylab="Interexceedance Times",cex=0.7)
  abline(v=qexp(1-Theta))
  abline(a = -qexp(1-Theta)/Theta, b=1/Theta)
}

declust <- function(x, ...)
{
  UseMethod("declust")
}

declust.extremalIndex <- function(Est)
{
  theCall <- match.call()
  
  C <- floor(Est$EstIntervals * Est$nExceed) + 1

  Times <- Est$interExceedTimes
  C <- min(C,length(Times)) # needed if short series and C < number of interexceedance times
  sTimes <- sort(Times, decr=TRUE)
  while(sTimes[C-1] == sTimes[C]){
    C <- C-1
  }

  clusters <- vector("list",sum(Times > sTimes[C])+1) # +1 since have one more exceedance than inter-exceedance times
  clusters[[1]] <- c(Est$exceedanceTimes[1], Est$thExceedances[1])
  Icluster <- 1

  for(i in 1:(Est$nExceed - 1)) {
    if(Times[i] > sTimes[C]){
      Icluster <- Icluster + 1
    }
    clusters[[Icluster]] <- rbind(clusters[[Icluster]],
                                  c(Est$exceedanceTimes[i+1], Est$thExceedances[i+1]))
  }
  
  clusterMaxima <- sapply(clusters, function(x) if(is.null(dim(x)))x[2] else max(x[,2]))
  
  res <- list(clusters = clusters,
              clusterMaxima = clusterMaxima,
              data = Est$data,
              threshold=Est$threshold,
              theCall=theCall)

  oldClass(res) <- "declustered"

  res
}

print.declustered <- function(x){
  print(x$theCall)
  cat("\nThreshold ",x$threshold,"\n")
  cat("Identified",length(x$clusters),"clusters.\n")
}

show.declustered <- print.declustered

plot.declustered <- function(x,ylab="Data",...){
  plot(x$data,xlab="",ylab=ylab)
  abline(h=x$threshold,col=2)
  for(i in 1:length(x$clusters)){
    clust <- x$clusters[[i]]
    if(is.null(dim(clust))){
      points(clust[1],clust[2],col=2)
    } else {
      points(clust,type="b",col=2)
    }
  }
}
