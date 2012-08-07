dgpd <- function(x, sigma, xi, u = 0, log.d=FALSE ){

    n <- length(x)

    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)
    u <- rep(u, length=n)

    .C(.c.dgpd, result=double(n), as.integer(n),
       as.double(x), as.double(sigma), as.double(xi),
       as.double(u), as.integer(log.d),
       NAOK=TRUE)$result
}

test.dgpd <- function(){
  evd.dgpd <- texmex:::.evd.dgpd

  myTest <- function(sig,xi,thresh,msg){
    myd <- sapply(1:nreps,function(i) dgpd(x[,i], sig[i], xi[i],u=thresh[i]))
    ed <- sapply(1:nreps, function(i) evd.dgpd(x[,i], loc=thresh[i], scale=sig[i], shape=xi[i]))
    checkEqualsNumeric(ed,myd,msg=msg)
    }

  set.seed(20101111)

#*************************************************************
# 6.12. Test dgpd. Note that .evd.dgpd is NOT vectorized.

  nreps <- 100
  nsim <- 1000
  p <- matrix(runif(2*nreps, -1, 1),ncol=2)
  p[, 1] <- p[, 1] + 1
  thresh <- rep(0,nreps)

  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))

  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="dgpd: random xi")

#*************************************************************
# 6.13. Test dgpd when some or all of xi == 0

  p[sample(1:nreps,nreps/2),2] <- 0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="dgpd: some zero xi")

  p[,2] <-  0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="dgpd: all zero xi")

#*************************************************************
# 6.14. Test vectorization of dgpd.

  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)
  thresh <- rnorm(nsim)

  x <- rgpd(nsim, sig, xi,u=thresh)
  myd <- dgpd(x, sig, xi,u=thresh)

  ed <- sapply(1:nsim, function(i) evd.dgpd(x[i], loc=thresh[i], scale=sig[i], shape=xi[i]))
  checkEqualsNumeric(ed,myd,msg="dgpd: vectorisation")

#*************************************************************
# 6.15 test log.d argument

  ld <- dgpd(x,sig,xi,u=thresh,log.d=TRUE)
  checkEqualsNumeric(myd,exp(ld),msg="dgpd: log density")
}
