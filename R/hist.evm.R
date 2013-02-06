hist.evmOpt <-
function(x, xlab, ylab, main, ...){
    a <- x$coefficients
#    a[1] <- exp(a[1])
    u <- x$threshold
    if (!is.finite(u)){ u <- min(x$data$y) }

    # FOLLOWING if BLOCK COMMENTED OUT TO ACCOUNT FOR DIFFERENCE
    # BETWEEN GEV AND GPD. MIGHT HAVE TO DO SOMETHING MORE
    # SENSIBLE LATER.
#    if(a[2] < 0){
#         UpperEndPoint <- u - a[1]/a[2]
#    }
#    else {
        UpperEndPoint <- Inf
#    }

    dat <- x$data$y
    dfun <- x$family$density

    h <- hist(dat, plot = FALSE)
    xx <- seq(u, min(UpperEndPoint, max(h$breaks)), length = 100)
    y <- dfun(xx, a, x)

    if (missing(xlab) || is.null(xlab)) xlab <- "Data"
    if (missing(ylab) || is.null(ylab)) ylab <- ""
    if (missing(main) || is.null(main)) main <- "Histogram and density"

    breaks <- seq(from=min(dat),to=max(dat),len=nclass.Sturges(dat)+1)

    hist(dat, prob = TRUE, ylim = c(0, max(y)),
         xlab=xlab, ylab=ylab, main=main, breaks = breaks, ...)
    lines(xx, y, col = 4)
    rug(dat)
    invisible()
}

