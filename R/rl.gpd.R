`rl.gpd` <-
function(object, alpha = .050,
         xlab, ylab, main ,
         pch= 1, col =2 , cex=.75, linecol = 4 ,
         cicol = 0, polycol = 15, smooth = TRUE ){

    # Define a helper function. Page 82 of Coles.
    gpd.delta <- function(a, m){
        # Jan's original code worked with sigma, not log(sigma)
        # Need to multiply the second component by exp(phi) = sigma.
        # This is not exact if a prior (penalty) function is used, but
        # the CI is approximate anyway.
        if (length(a) != 3){
            stop("covariates not allowed in return level plot")
        }
        
        out <- matrix(0, nrow=3, ncol=length(m))

        # Deal with exponential case
        if (a[3] == 0){
            out[1,] <- a[2] / a[1]
            out[2,] <- log(m * a[1])
        }
        else {
            out[1,] <- a[2] * m^a[3] * a[1]^(a[3] - 1)
            out[2,] <- 1/a[3] * ((m*a[1])^a[3] - 1) * exp(a[2]) # <--- mult by sigma = e^phi
            out[3,] <- -a[2] / (a[3]*a[3]) * ( (m * a[1] )^a[3] - 1 ) +
                       a[2] / a[3] * (m * a[1])^a[3] * log(m * a[1])
        } # Close else

        out

    } # Close gpd.delta

    a <- object$coefficients
    u <- object$threshold
    la <- object$rate
    n <- length(object$y) / la # Number of obs prior to thresholding

    mat <- object$cov
    xdat <- object$y
	
    a <- c(la, a)

    jj <- seq(-1, 3.75, by = 0.1)
    m <- unique( c(1/la, 10^jj) )

    q <- qgpd2(m, exp(a[2]), a[3], u, la)
    d <- t(gpd.delta(a = a, m = m))

    # Get covariance including P(over threshold) parameter
    mat <- matrix(c(la * (1 - la)/n, 0, 0,
                  0, mat[1,],
                  0, mat[2,]), ncol = 3)

    # Get (4.15) of Coles, page 82, adjusted for phi = log(sigma)
    v <- mahalanobis(d, center=c(0, 0, 0), cov=mat, inverted=TRUE)

    if (missing(xlab) || is.null(xlab)) { xlab <- "Return period" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Return level" }
    if (missing(main) || is.null(main)) { main <- "Return Level Plot" }

    plot(m, q,
         log = "x",
         type = "n",
         xlim = c(1, max(m)),
         ylim = c(u, max(xdat, q[q > u - 1] + qnorm(1-alpha/2) * sqrt(v)[q > u - 1])), 
         xlab = xlab, ylab = ylab, main = main)

    # Do polygon and CI lines
    if (smooth & length(xdat) > 2) {
        splo <- spline(m[q > u - 1],
                       q[q > u - 1] - qnorm(1-alpha/2) * sqrt(v)[q > u - 1] ,
                       200)
        sphi <- spline(m[q > u - 1],
                       q[q > u - 1] + qnorm(1-alpha/2) * sqrt(v)[q > u - 1] ,
                       200)
        if ( polycol != 0 ) {
            polygon( c( splo$x, rev( sphi$x ) ),
	             c( splo$y, rev( sphi$y ) ),
                     col = polycol ,
                     border=FALSE		    )
         } # Close if (polycol
         lines( splo$x, splo$y, col = cicol )
         lines( sphi$x, sphi$y, col = cicol )
    } # Close if (smooth &
    else{
        if (polycol != 0){
            polygon(c(m[q > u - 1], rev( m[q > u - 1])),
                    c(q[q > u - 1] - qnorm(1 - alpha/2) * sqrt(v)[q > u - 1],
                      rev(q[q > u - 1] + qnorm(1 - alpha/2) * sqrt(v)[q > u - 1])),
                    col=polycol,
                    border = FALSE) # Close polygon
        } # Close if
        else {
            lines(m[q > u - 1],
                  q[q > u - 1] + qnorm(1 - alpha/2) * sqrt(v)[q > u - 1],
                  col = cicol)
            lines(m[q > u - 1],
                  q[q > u - 1] - qnorm(1 - alpha/2) * sqrt(v)[q > u - 1],
                  col = cicol)
        }
    } # Close else
	
    lines(m[q > u - 1], q[q > u - 1],
          col = linecol[ 1 ] )

    # Add observed data to the plot
    ly <- length(object$y)
    points(1 / (1 - ((n - ly + 1):n) / (n + 1)), sort(xdat), pch=pch, col=col)
    box()
    invisible()
}

