pgpd <-
function(q, sigma, xi, u = 0, lower.tail=TRUE, log.p=FALSE ){
  ## first shift, rescale and cut at 0
  ## to handle q < u
  q <- pmax((q - u) / sigma, 0)

  ## now expand everything to the same length
  n <- max(length(q), length(xi))
  q <- rep(q, length.out=n)
  xi <- rep(xi, length.out=n)

  ## this handles negative xi properly
  xiq <- pmax(xi * q, -1)

  ## because we lack exprel and its inverse
  ## or, alternately, a lazy ifelse
  flag <- xi != 0
  expVal <- q
  expVal[flag] <- log1p(xiq[flag]) / xi[flag]
  expVal[is.na(xiq)] <- NA

  ## and now defer everything to pexp
  pexp(expVal, lower.tail=lower.tail, log.p=log.p)
}



test.pgpd <- function(){

  evd.pgpd <- texmex:::.evd.pgpd
  myTest <- function(sig,xi,thresh,msg){
    myp <- sapply(1:nreps,function(i) pgpd(x[,i], sig[i], xi[i],u=thresh[i]))
    ep <- sapply(1:nreps, function(i) evd.pgpd(x[,i], loc=thresh[i], scale=sig[i], shape=xi[i]))
    checkEqualsNumeric(ep,myp,msg=msg)
    }

  set.seed(20101111)

#*************************************************************
# 6.7. Test pgpd. Note that .evd.pgpd is NOT vectorized.

  nreps <- 100
  nsim <- 1000
  p <- matrix(runif(2*nreps, -1, 1),ncol=2)
  p[, 1] <- p[, 1] + 1
  thresh <- rep(0,nreps)

  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))

  myTest(sig=p[,1], xi=p[,2],thresh=thresh, msg="pgpd: random xi")

#*************************************************************
# 6.8. Test pgpd when some or all of xi == 0

  p[sample(1:nreps,nreps/2),2] <- 0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="pgpd: some zero xi")

  p[,2] <-  0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="pgpd: all zero xi")

#*************************************************************
# 6.9. Test vectorization of pgpd.

  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)
  thresh <- rnorm(nsim)

  x <- rgpd(nsim, sig, xi,u=thresh)
  myp <- pgpd(x, sig, xi,u=thresh)

  ep <- sapply(1:nsim, function(i)evd.pgpd(x[i], loc=thresh[i], scale=sig[i], shape=xi[i]))
  checkEqualsNumeric(ep,myp,msg="pgpd: vectorisation")

#*************************************************************
# 6.10 test log.p argument

  lp <- pgpd(x,sig,xi,u=thresh,log.p=TRUE)
  checkEqualsNumeric(myp,exp(lp),msg="pgpd: log probabilities")

#*************************************************************
# 6.11 test lower tail argument

  sp <- pgpd(x,sig,xi,u=thresh,lower.tail=FALSE)
  checkEqualsNumeric(myp,1-sp,msg="pgpd: lower tail")

  ## check pgpd when q < threshold
  upperProb <- pgpd(0, 1, 1, u=0.5, lower.tail=TRUE)
  checkEqualsNumeric(upperProb, 0, msg="pgpd: value below threshold (1)")

  lowerProb <- pgpd(0, 1, 1, u=0.5, lower.tail=FALSE)
  checkEqualsNumeric(upperProb, 0, msg="pgpd: value below threshold (2)")

  ## check pgpd when xi < 0 and value above upper limit

  xi <- -2.3
  upperProb <- pgpd(-2/xi, 1, xi, u=0, lower.tail=TRUE)
  checkEqualsNumeric(upperProb, 1, msg="pgpd: negative xi (1)")

  lowerProb <- pgpd(-2/xi, 1, xi, u=0, lower.tail=FALSE)
  checkEqualsNumeric(lowerProb, 0, msg="pgpd: negative xi (2)")
}

