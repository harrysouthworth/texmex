ppevm <-
function(object , nsim = 1000, alpha = .050,
		 xlab, ylab,  main, # labels and titles
		 pch=1, col = 2, cex = .75, linecol = 4 ,
		 intcol = 0, polycol=15){
    # Want parameters as a matrix with one row for passing
    # into family$rng and so on
    a <- t(object$coefficients)
    u <- object$threshold
    dat <- object$data$y

    pfun <- object$family$prob
    rfun <- object$family$rng

    if (missing(xlab) || is.null(xlab) ) { xlab <- "Model" }
    if (missing(ylab) || is.null(ylab) ) { ylab <- "Empirical" }
    if (missing(main) || is.null(main) ) { main <- "Probability Plot" }

    ModPoints <- ppoints(dat)

	# If doing the envelope, simulate, sort and get the quantiles
    if (nsim > 0){
        n <- length(dat)
        sim <- matrix(rfun(nsim * n, a, object), ncol = nsim)
        sim <- apply(sim, 2, sort)
        sim <- apply(sim, 2, pfun, param=a, model = object)
        sim <- apply(sim, 1, quantile, prob = c(alpha/2, 1 - alpha/2))
    }
    else { sim <- NULL }
  
    oldpar <- par(pty = "s"); on.exit(oldpar)
    plot(ModPoints, pfun(sort(dat), a, object),
         xlab = xlab, ylab = ylab, main = main, type = "n")

	# If doing the envelope, plot it before putting the data on 
    if (nsim > 0){
        if (polycol != 0)
            polygon(c(ModPoints, rev(ModPoints)), c(sim[1, ], rev(sim[2, ])),
                    col=polycol, border=FALSE)
        lines(ModPoints, sim[1, ], col = intcol)
        lines(ModPoints, sim[2, ], col = intcol)
    }

    abline(0, 1, col = linecol)
    points(ModPoints,
           pfun(sort(dat), a, object),
           pch = pch , col = col, cex = cex)
    box()
    invisible(sim)
}

