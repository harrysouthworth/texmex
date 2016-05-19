#' @method ggplot mrl
#' @export
ggplot.mrl <- function(data, xlab = "Threshold", ylab = "Mean excess", main=NULL,
                       fill="light blue", col="blue",
                       addNexcesses=TRUE, textsize=4, ...){
  x <- data
  data <- x$data
  x <- x$mrl
  d <- data.frame(th = x[, "threshold"],
                  mrl = x[, "MRL"],
                  xl = x[, "lo"],
                  xu = x[, "hi"])

  k <- !is.na(d$xl)
  poly <- data.frame(x=c(d$th, rev(d$th)), y=c(d$xl, rev(d$xu)))
  poly <- poly[c(k, rev(k)), ]

  p <- ggplot(poly, aes(x, y)) +
    geom_polygon(fill=fill, alpha=.5) +
    geom_line(data=d, aes(th, mrl), color=col) +
    scale_x_continuous(xlab) +
    scale_y_continuous(ylab) +
    ggtitle(main)

  if (addNexcesses)
    p <- addExcesses(p, poly$x, poly$y, data=data, u=d$th,
                     textsize=textsize)

  p
}

#' @method ggplot gpdRangeFit
#' @export
ggplot.gpdRangeFit <- function(data, xlab = "Threshold", ylab = NULL, main = NULL,
                               fill="orange", col="blue",
                               addNexcesses = TRUE, textsize=4, ...){
  if (missing(ylab)) {
    ylab <- c(expression(hat(phi)[m]), expression(hat(xi)))
  }  else if (length(ylab) != 2) {
    stop("length of ylab should be 2")
  }
  if (!missing(main) && length(main) != 2) {
    stop("length of main should be 2")
  }

  x <- data
  data <- data$data

  p <- vector("list", 2)

  for (i in 1:2) {
    #        yl <- range(x$hi[, i], x$lo[, i])

    d <- data.frame(th=x$th, par=x$par[, i])
    poly <- data.frame(x=c(x$th, rev(x$th)), y=c(x$lo[, i], rev(x$hi[, i])))

    p[[i]] <- ggplot(poly, aes(x, y)) +
      geom_polygon(fill=fill, alpha=.5) +
      geom_line(data=d, aes(th, par), color=col) +
      scale_x_continuous(xlab) +
      scale_y_continuous(ylab[i]) +
      theme(axis.title.y=element_text(angle=0)) +
      if (!missing(main)) ggtitle(main[i])

    if (addNexcesses)
      p[[i]] <- addExcesses(p[[i]], poly$x, poly$y, data=data, u=u, textsize=textsize)
  } # Close for
  p
}
