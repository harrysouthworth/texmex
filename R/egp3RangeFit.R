egp3RangeFit <-
function (data, umin=quantile(data, .05), umax=quantile(data, .95),
          nint = 10,
          penalty="gaussian", priorParameters=NULL, alpha=.05) {

  if (umin < 0)
    stop("umin < 0: data must be non-negative. Add a constant or increase umin and try again")

  m <- s <- hi <- lo <- rep(0, nint)
  u <- seq(umin, umax, length = nint)
  qz <- qnorm(1 - alpha/2)

  for (i in 1:nint) {
    z <- evm(data, th=u[i], penalty=penalty, priorParameters=priorParameters, family=egp3)
    m[i] <- z$coefficients[3]
    s[i] <- z$se[3]
  }

  # egp3 family works with labmda = log(kappa)
  hi <- exp(m + qz * s)
  lo <- exp(m - qz * s)

  res <- list(th=u, par=exp(m) , hi=hi, lo=lo, data=data)
  oldClass(res) <- 'egp3RangeFit'
  res
}

print.egp3RangeFit <- function(x, ...){
  print(cbind(threshold=x$th, kappa=x$par, lo=x$lo, hi=x$hi))
  invisible()
}

plot.egp3RangeFit <- function(x, xlab="Threshold", ylab="kappa",
                              main=NULL, addNexcesses=TRUE, ...){

  yl <- range(x$hi, x$lo, 1)
  plot(x$th, x$par, ylim = yl, type = "b",
       xlab=xlab, ylab=ylab, main=main, log="y", ...)
  for (j in 1:length(x$th)){
    lines(c(x$th[j], x$th[j]), c(x$hi[j], x$lo[j]))
  }
  if(addNexcesses){
    axis(3, at=axTicks(1), labels=sapply(axTicks(1), function(u) sum(x$data > u)), cex=0.5)
    mtext("# threshold excesses")
  }
  abline(h=1, lty=2)
  
  invisible()
}