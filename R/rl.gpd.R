`rl.gpd` <-
function(object, alpha = .050, RetPeriodRange=NULL,
         xlab, ylab, main,
         pch= 1, col =2 , cex=.75, linecol = 4 ,
         cicol = 0, polycol = 15, smooth = TRUE ){

    # Define a helper function. Page 82 of Coles, note that our parameters are phi=log sigma and xi so implementation slightly different.
    gpd.delta <- function(a, m){
        # This is not exact if a prior (penalty) function is used, but
        # the CI is approximate anyway.
        if (length(a) != 3){
            stop("covariates not allowed in return level plot")
        }
        
        out <- matrix(0, nrow=3, ncol=length(m))
        
        if (a[3] == 0){ # exponential case
            out[1,] <- exp(a[2]) / a[1]
            out[2,] <- exp(a[2]) * log(m * a[1])
        } else {
            out[1,] <- exp(a[2]) * m^a[3] * a[1]^(a[3] - 1)
            out[2,] <- exp(a[2]) / a[3] * ((m*a[1])^a[3] - 1) 
            out[3,] <- -exp(a[2]) / (a[3]*a[3]) * ( (m * a[1] )^a[3] - 1 ) +
                       exp(a[2]) / a[3] * (m * a[1])^a[3] * log(m * a[1])
        } 

        out
    } 

    a <- object$coefficients
    u <- object$threshold
    la <- object$rate # rate of threshold excess
    n <- length(object$y) / la # Number of obs prior to thresholding

    xdat <- object$y
	
    a <- c(la, a)

    if( is.null(RetPeriodRange) ){
      jj <- seq(-1, max(3.75,log10(n)),by=0.1)
    }  else {
      jj <- seq(log10(RetPeriodRange[1]),log10(RetPeriodRange[2]),length=400)
    }
    m <- unique( c(1/la, 10^jj) )

    xm <- qgpd2(m, exp(a[2]), a[3], u, la)
    dxm <- t(gpd.delta(a = a, m = m))

    # Get covariance including P(over threshold) parameter
    V <- matrix(c(la * (1 - la)/n, 0, 0,
                  0, object$cov[1,],
                  0, object$cov[2,]), ncol = 3)

    # Get (4.15) of Coles, page 82, adjusted for phi = log(sigma)
    vxm <- mahalanobis(dxm, center=c(0, 0, 0), cov=V, inverted=TRUE)

    if (missing(xlab) || is.null(xlab)) { xlab <- "Return period" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Return level" }
    if (missing(main) || is.null(main)) { main <- "Return Level Plot" }

    plot(m, xm,
         log = "x",
         type = "n",
         xlim=range(m),
         ylim=range(c(xdat, xm[xm > u - 1] + qnorm(1-alpha/2) * sqrt(vxm)[xm > u - 1])), 
         xlab = xlab, ylab = ylab, main = main)

    # Do polygon and CI lines
    U <- u - abs(u/100)
    if (smooth & length(xdat) > 2) {
        splo <- spline(log(m[xm > U]),
                       xm[xm > U] - qnorm(1-alpha/2) * sqrt(vxm)[xm > U] ,
                       200)
        sphi <- spline(log(m[xm > U]),
                       xm[xm > U] + qnorm(1-alpha/2) * sqrt(vxm)[xm > U] ,
                       200)
        if ( polycol != 0 ) {
            polygon( exp(c( splo$x, rev( sphi$x ) )),
	             c( splo$y, rev( sphi$y ) ),
                     col = polycol ,
                     border=FALSE		    )
         } # Close if (polycol
         lines( exp(splo$x), splo$y, col = cicol )
         lines( exp(sphi$x), sphi$y, col = cicol )
    } else{
        if (polycol != 0){
            polygon(c(m[xm > U], rev( m[xm > U])),
                    c(xm[xm > U] - qnorm(1 - alpha/2) * sqrt(vxm)[xm > U],
                      rev(xm[xm > U] + qnorm(1 - alpha/2) * sqrt(vxm)[xm > U])),
                    col=polycol,
                    border = FALSE) # Close polygon
        } else {
            lines(m[xm > U],
                  xm[xm > U] + qnorm(1 - alpha/2) * sqrt(vxm)[xm > U],
                  col = cicol)
            lines(m[xm > U],
                  xm[xm > U] - qnorm(1 - alpha/2) * sqrt(vxm)[xm > U],
                  col = cicol)
        }
    } 
	
    lines(m[xm > U], xm[xm > U], col = linecol[ 1 ] )

    # Add observed data to the plot
    ly <- length(xdat)
    points(1 / (1 - ((n - ly + 1):n) / (n + 1)), sort(xdat), pch=pch, col=col)
    box()
    invisible()
}

