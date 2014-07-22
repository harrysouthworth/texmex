egp3RangeFit <-
function (data, umin=quantile(data, .05), umax=quantile(data, .95),
          nint = 10,
          penalty="gaussian", priorParameters=NULL, alpha=.05) {
    
  m <- s <- hi <- lo <- rep(0, nint)
  u <- seq(umin, umax, length = nint)
  qz <- qnorm(1 - alpha/2)

  for (i in 1:nint) {
    z <- evm(data, th=u[i], penalty=penalty, priorParameters=priorParameters, family=egp3)
    m[i] <- z$coefficients[3]
    s[i] <- sqrt(z$se[3])
  }

  # egp3 family works with labmda = log(kappa)
  hi[i, ] <- exp(m[i] + qz * s)
  lo[i, ] <- exp(m[i] - qz * s)

  res <- list(th=u, par=m , hi=hi, lo=lo, data=data)
  oldClass(res) <- 'egp3RangeFit'
  res
}

print.egp3RangeFit <- function(x, ...){
  cbind(threshold=x$u, kappa=x$m, lo=x$lo, hi=x$hi)
  invisible()
}

