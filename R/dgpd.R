`dgpd` <-
function(x, sigma, xi, u = 0, log.d=FALSE ){

    n <- length(x)

    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)
    u <- rep(u, length=n)

    if (all(xi == 0)){
      res <- dexp(x-u, 1/sigma, log=TRUE)
    }
    else if (any(xi == 0)){
      res <- numeric(n)
	    wh <- xi == 0
	    res[wh] <- dexp(x[wh] - u[wh], 1/sigma[wh], log=TRUE)
	    res[!wh] <- logb(1 + (xi[!wh] * (x[!wh] - u[!wh]))/sigma[!wh] ) * ( -1/xi[!wh] - 1) - log(sigma[!wh])
    }
    else{
      res <- logb(1 + (xi * (x - u))/sigma ) * ( -1/xi - 1) - log(sigma)
    }

    res <- ifelse( exp(res) >=0, res, 0 )
	
    if (!log.d){
        res <- exp(res)
    }
    res
}

test.dgpd <- function(){

  myTest <- function(sig,xi,thresh,msg){
    myd <- sapply(1:nreps,function(i) dgpd(x[,i], sig[i], xi[i],u=thresh[i]))
    ed <- sapply(1:nreps, function(i) .evd.dgpd(x[,i], loc=thresh[i], scale=sig[i], shape=xi[i]))
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

  x <- sapply(1:nreps,function(i)rgpd(nsim,sig=p[i,1],xi=p[i,2],u=thresh[i]))

  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="dgpd: random xi")

#*************************************************************
# 6.13. Test dgpd when some or all of xi == 0

  p[sample(1:nreps,nreps/2),2] <- 0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sig=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="dgpd: some zero xi")

  p[,2] <-  0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sig=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="dgpd: all zero xi")

#*************************************************************
# 6.14. Test vectorization of dgpd.

  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)
  thresh <- rnorm(nsim)
  
  x <- rgpd(nsim, sig, xi,u=thresh)
  myd <- dgpd(x, sig, xi,u=thresh)
  
  ed <- sapply(1:nsim, function(i).evd.dgpd(x[i], loc=thresh[i], scale=sig[i], shape=xi[i]))
  checkEqualsNumeric(ed,myd,msg="dgpd: vectorisation")

#*************************************************************
# 6.15 test log.d argument
  
  ld <- dgpd(x,sig,xi,u=thresh,log.d=TRUE)
  checkEqualsNumeric(myd,exp(ld),msg="dgpd: log density")
}
