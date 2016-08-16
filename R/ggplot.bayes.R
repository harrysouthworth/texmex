ggdensplots <- function(x, fill="dark blue", col="light blue"){
    v <- names(coef(x))
    n <- length(v)
    p <- vector("list", length=n)

    for (i in 1:n){
        d <- data.frame(x$param[, i])
        names(d) <- "x"
        p[[i]] <- ggplot(data=d, aes(x=x)) +
                     stat_density(fill=fill,alpha=0.5) +
                     scale_x_continuous(v[i]) +
                     scale_y_continuous("")
    }
    p
}

ggacfplots <- function(x, fill="orange"){
    v <- names(coef(x))
    n <- length(v)
    p <- vector("list", length=n)
    thin <- x$thin

    for (i in 1:n){
        acz <- acf(x$param[, i], plot=FALSE)
        acd <- data.frame(lag=acz$lag, acf=acz$acf, xend=acz$lag, yend=rep(0, length(acz$lag)))

        p[[i]] <- ggplot(acd, aes(lag, acf)) +
                      geom_area(fill=fill) +
                      geom_segment(color="dark blue", aes(x=lag, y=acf, xend=xend, yend=yend), alpha=0.5) +
                      scale_x_continuous("Lag") +
                      scale_y_continuous("ACF") +
                      ggtitle(paste0("ACF for ", v[i], "\n(thin = ", thin, ")"))
    }
    p
}

ggtraceplots <- function(x, trace="light blue", mean="blue", burn="orange"){
    v <- names(coef(x))
    n <- length(v)
    p <- vector("list", length=n)

    for (i in 1:n){
        cm <- cumsum(x$chains[, i]) / 1:nrow(x$chains)
        d <- data.frame(x=1:nrow(x$chains), cm=cm, p=x$chains[, i])
        rects <- data.frame(x=c(0, x$burn), xend=c(x$burn, nrow(x$chains)),
                            col=as.character(c(burn, "light grey")))

        p[[i]] <- ggplot() +
                      geom_line(data=d, aes(x, p), color=trace) +
                      geom_line(data=d, aes(x, cm), color=mean, size=1.5) +
                      geom_rect(data=rects, aes(xmin=x, xmax=xend,
                                                ymin=-Inf, ymax=Inf),
                                                fill=c(burn, "grey90"), alpha=.2) +
                      scale_x_continuous("Step number") +
                      scale_y_continuous(paste(v[i], "\n& cumulative mean"))

    }
    p
}

#' Diagnostic plots for the Markov chains in an evmSim object
#' @param data An object of class 'evmSim'.
#' @param which.plots Which plots to produce. Density plots correspond
#'     to 1, trace plots of the Markov chains to 2 and autocorrelation
#'     function plots to 3.
#' @param denscol Colour for the density plots. Defaults to 'dark blue'.
#' @param acfcol Colour for the ACF plots. Defaults to 'light blue'.
#' @param plot.it Whether or not to actually print the plots. Defaults
#'     to \code{plot.it=TRUE}.  If \code{plot.it=FALSE}, you might
#'     want to control the layout. Do this with
#'     \code{do.call("grid.arrange", c(plots, ncol=2))}, for example,
#'     where \code{plots} is the objected returned by
#'     \code{ggplot.evmSim}.
#' @param mapping,environment ignored
#' @param ... Additional arguments to \code{ggplot}, currently unused.
#' @aliases ggtraceplots ggdensplots ggacfplots
#' @keywords hplot
#' @method ggplot evmSim
#' @export
ggplot.evmSim <- function(data=NULL, mapping, which.plots=1:3, denscol="dark blue", acfcol="light blue", plot.it=TRUE,
                          ..., environment){
    d <- if (1 %in% which.plots) ggdensplots(data, fill=denscol)
         else NULL
    tr <- if (2 %in% which.plots) ggtraceplots(data)
          else NULL
    a <- if (3 %in% which.plots) ggacfplots(data, fill=acfcol)
         else NULL
    res <- c(d, tr, a)
    if (plot.it) do.call("grid.arrange", c(res, ncol=ncol(data$param)))
    invisible(res)
}
