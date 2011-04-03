#---------------------------------------------------------------------------
#Function MCS() evaluates the MCS wrt to Schmid et al. 2006 equation (17)
#First and second args are
#X (matrix) and p (vector of probabilities).
#Matrix X has to be i*j where i=1,...,n, and j=1,...,d and d,n
#denote dimension and length, respectively.	
#---------------------------------------------------------------------------

.MCSlower <- function(U,p)
  {
    d <- dim(U)[1]
    n <- dim(U)[2]
    res1         <- p-U
    res1[res1<0] <- 0
    res2         <- apply(res1,2,prod)
    res3         <- sum(res2)
    res4         <- ( (1/n)*res3-((p^2)/2)^(d) )/
      ( (p^(d+1))/(d+1) -((p^2)/2)^(d))
    return(res4)
  }

MCS <- function(X,p=seq(.1, .9, by=.1)) {
    theCall <- match.call()
     
    X <- t(X) # Yiannis's original code had variables as rows
    U <- t(apply(X,1,edf)) #transpose cause apply transposes g(X), g:edf
    n    <- length(p)
	res <- sapply(p, .MCSlower, U=U)

    res <- list(mcs=res, p=p, call=theCall)
    oldClass(res) <- "MCS"
    res
  }

plot.MCS <- function(x, xlab="p", ylab= "MCS", ...){
   plot(x$p, x$mcs, type="l", xlab=xlab, ylab=ylab, ...)
   invisible()
}

print.MCS <- function(x, ...){
    print(x$call)
    cat("Multivariate conditional Spearman's rho.\n\n", sep = "")
    res <- x$mcs
    names(res) <- x$p
    print(res)
    invisible(res)
}
show.MCS <- print.MCS

summary.MCS <-  function(object, ...){
    print(object$call)
    cat("Multivariate conditional Spearman's rho.\n\n", sep = "")
    res <- object$mcs
    names(res) <- object$p
    print(res)
    invisible(res)
}

#------------------------------------------------
#Bootstrap
#------------------------------------------------

bootMCS <- function(X,p=seq(.1, .9, by=.1),R=100, trace=10) {
   theCall <- match.call()
   bfun <- function(i, data, p, trace){
       if (i %% trace == 0){ cat("Replicate", i, "\n") }
       d <- data[sample(1:nrow(data), replace=TRUE),]
       MCS(d, p)$mcs
   }

   res <- sapply(1:R, bfun, data=X, p=p, trace=trace)
   res <- list(replicates=res, p=p, R=R, call=theCall)
   oldClass(res) <- "bootMCS"
   invisible(res)
}

plot.bootMCS <- function(x, xlab="p", ylab= "MCS",alpha=.05, ylim, ...){
   m <- rowMeans(x$replicates)
   ci <- apply(x$replicates, 1, quantile, prob=c(1-alpha/2, alpha/2))
   if (missing(ylim)){ ylim <- range(ci) }
   plot(x$p, m, type="l", ylim=ylim,
        xlab=xlab, ylab=ylab,
        sub=paste(100*(1-alpha), "% interval. ", x$R, " bootstrap samples were performed", sep=""))
   lines(x$p, ci[1,], lty=2)
   lines(x$p, ci[2,], lty=2)
   invisible(ci)
}

print.bootMCS <- function(x, ...){
    print(x$call)
    cat("Multivariate conditional Spearman's rho.\n", x$R, " bootstrap samples were performed.\n\n",
        sep = "")
     
    m <- rowMeans(x$replicates)
    names(m) <- x$p
    print(m)
    invisible(m)
}
show.bootMCS <- print.bootMCS

summary.bootMCS <- function(object, alpha=.05, ...){
    cat("Multivariate conditional Spearman's rho.\n", object$R, " bootstrap samples were performed.\n\n",
        sep = "")
    m <- rowMeans(object$replicates)
    ci <- apply(object$replicates, 1, quantile, prob=c(alpha/2, 1 - alpha/2))
    res <- cbind(m, t(ci))
    dimnames(res) <- list(object$p, c("Mean", rownames(ci)))
    res
}

test.MCS <- function(){
  myMCS <- function(x,p){
# First and second args are
# x (dxn matrix) and p (vector of probabilities).

    n <- dim(x)[2]
    d <- dim(x)[1]
    u <- t(apply(x,1,rank))/(n+1)# schmid and schmidt use n not n+1
    
    rho <- numeric(length(p))
    for(k in 1:length(p)) {
      Diff <- p[k] - u
      Diff[Diff<0] <- 0
      Prod <- apply(Diff,2,prod)
      
      num <- mean(Prod) - ((p[k]^2)/2)^d
      den <- (p[k]^(d+1)) / (d+1) - (0.5*p[k]^2)^d
      rho[k] <- num/den
    }

    rho
   }
# simulated data - dimension 2
   n <- 1000
   by <- 0.01
   p <- seq(by,1-by,by=by)
   data <- rbind(rnorm(n),rnorm(n))

   tmRl <- MCS(t(data),p)
   myRl <- myMCS(data,tmRl$p)
   
   checkEqualsNumeric(myRl,tmRl$mcs,msg="MCS: independent normal data")
   checkEqualsNumeric(p,tmRl$p, msg="MCS: mathching p argument")
   
# winter air pollution data - dimension 5
   tmWinterMCS <- MCS(winter,p)
   myWinterMCS <- myMCS(t(winter),p)
   checkEqualsNumeric(myWinterMCS,tmWinterMCS$mcs,msg="MCS: winter air pollution data")
   
# summer airpollution data - dimension 5
   tmSummerMCS <- MCS(summer,p)
   mySummerMCS <- myMCS(t(summer),p)
   checkEqualsNumeric(mySummerMCS,tmSummerMCS$mcs,msg="MCS: summer air pollution data")
}
