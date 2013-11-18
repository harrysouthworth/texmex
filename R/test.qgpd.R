test.qgpd <-
function(){
  
  set.seed(201110101)
  evd.qgpd <- sombrero:::.evd.qgpd
  myTest <- function(sig,xi,thresh,msg){
    myq <- sapply(1:nreps,function(i) qgpd(x[,i], sig[i], xi[i], u=thresh[i]))
    myp <- sapply(1:nreps,function(i) pgpd(myq[,i], sig[i], xi[i], u=thresh[i]))
    eq <- sapply(1:nreps, function(i) evd.qgpd(x[,i], loc=thresh[i], scale=sig[i], shape=xi[i]))
    checkEqualsNumeric(eq,myq,msg=paste(msg,"test using .evd.qgpd"))
    checkEqualsNumeric(x,myp,msg=paste(msg,"test using qgpd"))
  }
  
  #*************************************************************
  # 6.4.0 Test exception for out of range probabilties
  op <- options()
  options(show.error.messages=FALSE)
  checkException(qgpd(1.5,1,0,2),msg="qgpd: exception for out of range prob")
  checkException(qgpd(-1,1,0,2),msg="qgpd: exception for out of range prob")
  options(op)
  
  #*************************************************************
  # 6.4. Test qgpd. Note that .evd.qgpd is NOT vectorized.
  
  nreps <- 100
  nsim <- 1000
  p <- matrix(runif(2*nreps, -1, 1),ncol=2) 
  p[, 1] <- p[, 1] + 1
  thresh <- rep(0,nreps) 
  x <- matrix(runif(nreps*nsim), nrow=nsim)
  
  myTest(sig=p[,1], xi=p[,2],thresh=thresh,msg="qgpd: random xi")
  
  #*************************************************************
  # 6.5. Test qgpd when some or all of xi == 0. Note that .evd.qgpd is NOT vectorized.
  
  p[sample(1:nreps,nreps/2),2] <- 0
  myTest(sig=p[,1], xi = p[,2], thresh=thresh,msg="qgpd: some zero xi")
  p[,2] <-  0
  myTest(sig=p[,1], xi = p[,2], thresh=thresh,msg="qgpd: all zero xi")
  
  #*************************************************************
  # 6.6. Test vectorization of qgpd. Note that .evd.qgpd is NOT vectorized.
  
  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)
  thresh <- rnorm(nsim)
  
  x <- runif(nsim)
  
  myq <- qgpd(x, sig, xi, thresh)
  eq <- sapply(1:nsim, function(i)evd.qgpd(x[i], loc=thresh[i], scale=sig[i], shape=xi[i]))
  
  checkEqualsNumeric(eq,myq,msg="qgpd: vectorisation")
  
  #*************************************************************
  # 6.6a Test log.p argument
  
  lq <- qgpd(log(x), sig,xi,thresh,log.p=TRUE)
  
  checkEqualsNumeric(myq, lq, msg="qgpd: log.p=TRUE")
  
  #*************************************************************
  # 6.6a Test log.p argument
  
  LTq <- qgpd(1-x, sig,xi,thresh, lower.tail=FALSE)
  
  checkEqualsNumeric(myq, LTq, msg="qgpd: lower.tail=FALSE")
  
}
