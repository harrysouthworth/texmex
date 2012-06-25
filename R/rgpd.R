rgpd <- function(n, sigma, xi, u = 0) {
    # Check parameter vectors have correct length
    sigma <- as.double(rep(sigma, length=n))
    xi <- as.double(rep(xi, length=n))
    u <- as.double(rep(u, length=n))
    result <- double(n)
    # and pass to the underlying C
    .C(.c.rgpd, result=result, as.integer(n), sigma, xi, u)$result
  }

test.rgpd <- function(){
  evd.rgpd <- texmex:::.evd.rgpd
  myTest <- function(seed,p, thresh,msg=""){
      set.seed(seed)
      x <- sapply(1:nreps, function(i)rgpd(nsim, p[i,1], p[i,2], u=thresh[i]))
      set.seed(seed)  
      ex <- sapply(1:nreps, function(i)evd.rgpd(nsim, loc=thresh[i], scale=p[i,1],shape=p[i,2]))
      checkEqualsNumeric(ex,x,msg=msg)
      }
  seed <- 20101111
#*************************************************************
# 6.1. Test rgpd. Note that .evd.rgpd is NOT vectorized.

  nreps <- 100
  nsim <- 1000
  p <- matrix(runif(2*nreps, -1, 1),ncol=2) 
  p[, 1] <- p[, 1] + 1
  thresh <- rep(0,nreps) 
  
  myTest(seed,p,thresh=thresh,msg="rgpd: random parameters, zero threshold")
  
#*************************************************************
# 6.1a Test rgpd with non-zero threshold. Note that .evd.rgpd is NOT vectorized.

  nonZeroThresh <- rnorm(nreps)
  myTest(seed,p,thresh=nonZeroThresh,msg="rgpd: Non-zero threshold")
  
#*************************************************************
# 6.2. Test rgpd when some or all xi == 0. Note that .evd.rgpd is NOT vectorized.

  p[sample(1:nreps,nreps/2),2] <- 0
  myTest(seed,p,thresh=thresh,msg="rgpd: some zero xi")
  p[,2] <- 0
  myTest(seed,p,thresh=thresh,msg="rgpd: all zero xi")

#*************************************************************
# 6.3. Test vectorization of rgpd. .evd.rgpd is NOT vectorized

  sig <- runif(nreps, 0, 2)
  xi <- runif(nreps)

  set.seed(seed)
  x <- rgpd(nreps, sig, xi)
  set.seed(seed)
  ex <- sapply(1:nreps, function(i)evd.rgpd(1, loc=0, scale=sig[i], shape=xi[i]))
  
  checkEqualsNumeric(ex, x, msg="rgpd: vectorisation")
}

