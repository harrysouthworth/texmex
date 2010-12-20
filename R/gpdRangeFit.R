`gpdRangeFit` <-
function (data, umin=quantile(data, .05), umax=quantile(data, .95),
          nint = 10, 
          penalty="gaussian", priorParameters=NULL, alpha=.05,
          xlab="Threshold", ylab=NULL,
		  main=NULL, ... ) {
	if ( missing( ylab ) ){
		ylab = c( "log(scale)", "shape" )
	}
	else if ( length( ylab ) != 2 ){
		stop( "length of ylab should be 2" )
	}

	if ( !missing( main ) && length( main ) != 2 ){
		stop( "length of main should be 2" )
	}

    m <- s <- up <- ul <- matrix(0, nrow = nint, ncol = 2)
    u <- seq(umin, umax, length = nint)
    qz <- qnorm(1-alpha/2)
    for (i in 1:nint) {
        z <- gpd(data, th=u[i], penalty=penalty, priorParameters=priorParameters)
        m[i, ] <- z$coefficients
        m[i, 1] <- m[i, 1] - m[i, 2] * u[i]
        d <- matrix(c(1, -u[i]), ncol = 1)
        v <- t(d) %*% z$cov %*% d
        s[i, ] <- sqrt( diag( z$cov ) )
        s[i, 1] <- sqrt(v)
        
        up[i, ] <- m[i, ] + qz * s[i, ]
        ul[i, ] <- m[i, ] - qz * s[i, ]
    }
    names <- c("Modified Scale", "Shape")
    for (i in 1:2) {
        um <- max(up[, i])
        ud <- min(ul[, i])
        plot(u, m[, i], ylim = c(ud, um), type = "b",
			xlab=xlab, ylab=ylab[i], main=main[i], ...)
        for (j in 1:nint) lines(c(u[j], u[j]), c(ul[j, i], up[j, i]))
    }
    invisible()
}

test.gpdRangeFit <- function(){
  par(mfrow=c(2,1))
  res <- gpdRangeFit(rain, umin=0, umax=50, nint=20, pch=16, main=c("Figure 4.2 of Coles (2001)",""))
  checkEquals(res,NULL,msg="gpdRangeFit: check execution")
}
