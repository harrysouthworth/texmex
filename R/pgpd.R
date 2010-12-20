pgpd <-
function(q, sigma, xi, u = 0, lower.tail=TRUE, log.p=FALSE ){

    q <- (q - u) / sigma
    n <- length(q)

    xi <- rep(xi, length=n)
    sigma <- rep(sigma, length=n)
    u <- rep(u, length=n)

    if (all(xi == 0)){
        res <- pexp(q, log.p=TRUE, lower.tail = FALSE)
    }
    else if (any(xi == 0)){
        res <- numeric(n)
        wh <- xi == 0
        res[wh] <- pexp(q[wh], log.p=TRUE , lower.tail=FALSE)
        res[!wh] <- log(1 + xi[!wh]*q[!wh]) * (-1/xi[!wh])
    }
    else {
        res <- log(1 + xi * q) * (-1/xi)
    }

    if (!log.p){
        res <- exp(res) # survivor function
        if (lower.tail){
            res <- 1 - res
        }
    } # Close if (!log.p
    else { # if want log
        if(lower.tail){
            res <- log(1 - exp(res))
        }
    }

	res
}

test.pgpd <- function(){

  myTest <- function(sig,xi,thresh,msg){
    myp <- sapply(1:nreps,function(i) pgpd(x[,i], sig[i], xi[i],u=thresh[i]))
    ep <- sapply(1:nreps, function(i) .evd.pgpd(x[,i], loc=thresh[i], scale=sig[i], shape=xi[i]))
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

  x <- sapply(1:nreps,function(i)rgpd(nsim,sig=p[i,1],xi=p[i,2],u=thresh[i]))

  myTest(sig=p[,1], xi=p[,2],thresh=thresh, msg="pgpd: random xi")

#*************************************************************
# 6.8. Test pgpd when some or all of xi == 0

  p[sample(1:nreps,nreps/2),2] <- 0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sig=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="pgpd: some zero xi")

  p[,2] <-  0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sig=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="pgpd: all zero xi")

#*************************************************************
# 6.9. Test vectorization of pgpd.

  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)
  thresh <- rnorm(nsim)
  
  x <- rgpd(nsim, sig, xi,u=thresh)
  myp <- pgpd(x, sig, xi,u=thresh)
  
  ep <- sapply(1:nsim, function(i).evd.pgpd(x[i], loc=thresh[i], scale=sig[i], shape=xi[i]))
  checkEqualsNumeric(ep,myp,msg="pgpd: vectorisation")

#*************************************************************
# 6.10 test log.p argument
  
  lp <- pgpd(x,sig,xi,u=thresh,log.p=TRUE)
  checkEqualsNumeric(myp,exp(lp),msg="pgpd: log probabilities")
  
#*************************************************************
# 6.11 test lower tail argument

  sp <- pgpd(x,sig,xi,u=thresh,lower.tail=FALSE)
  checkEqualsNumeric(myp,1-sp,msg="pgpd: lower tail")
}

