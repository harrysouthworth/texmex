mexMonteCarlo <- function(nSample,mexList,mult=10){
#set up
  d <- length(mexList)
  data <- mexList[[1]]$margins$data
  margins <- mexList[[1]]$dependence$margins
  nData <- dim(data)[1]
#Generate our Monte Carlo sample from the original dataset. 
  which <- sample(1:nData, size=nSample, replace = TRUE)
  MCsampleOriginal <- data[which,]
#Transform the original dataset to the Laplace scale. 
  dataLaplace <- mexTransform(mexList[[1]]$margins, margins = margins, method = "mixture")$transformed
  MCsampleLaplace <- dataLaplace[which,]
#Identify maximum components. 
  whichMax <- apply(MCsampleLaplace,1,which.max)
# Now identify which of the maximal components lie above their associated conditional dependence model thresholds.
  dth <- sapply(mexList,function(l)l$dependence$dth)
  dqu <- sapply(mexList,function(l)l$dependence$dqu)
  whichMaxAboveThresh <- sapply(1:nSample,function(i)MCsampleLaplace[i,whichMax[i]] >= dth[whichMax[i]])
# generate large samples from each of the Conditional models,
  mexKeep <- lapply(1:d,function(i){
      mc <- predict.mex(mexList[[i]],pqu=dqu[i],nsim=nSample*d*mult)
      mc$data$simulated[mc$data$CondLargest,c(i,c(1:d)[-i])]})
# Replace original sample by samples from conditional models. 
  nR <- rep(0,d)
  names(nR) <- names(data)
  for(i in 1:d){
      replace <- whichMax == i & whichMaxAboveThresh
      nReplace <- sum(replace)
      if(nReplace > 0){
         nR[i] <- nReplace
         MCsampleOriginal[replace,] <- as.matrix(mexKeep[[i]])[1:nReplace,]
      }
  }

  list(nR=nR, MCsample=MCsampleOriginal, whichMax=whichMax, whichMaxAboveThresh=whichMaxAboveThresh)
}

mexAll <- function(data,mqu,dqu) {
    d <- dim(data)[2]
    res <- lapply(1:d,function(i){mex(data,which=i,mqu=mqu,dqu=dqu[i])})
    names(res) <- names(data)
    oldClass(res) <- "mexList"
    res
}

print.mexList <- function(x, ...){
    print(x[[1]])
    for(i in 2:length(x)){
        cat("\n______\n")
        print(x[[i]][[2]])
    }
}
