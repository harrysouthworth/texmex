`cmvxDependence` <-
function (x, which, gth, gqu)
{
   theCall <- match.call()
   if (class(x) != "migpd")
       stop("you need to use an object created by migpd and passed through cmvxGumbel")
   else if (is.null(x$gumbel))
       stop("you need to pass the object through cmvxGumbel")
   if (missing(which)) {
       cat("Missing 'which'. Conditioning on", dimnames(x$gumbel)[[2]][1],
           "\n")
       which <- 1
   }
   else if (length(which) > 1)
       stop("which must be of length 1")
   else if (is.character(which))
       which <- match(which, dimnames(x$gumbel)[[2]])
   if (missing(gth) & missing(gqu)) {
       cat("Assuming same quantile for thesholding as was used on raw data...\n")
       gqu <- x$qu[which]
   }
   else if (missing(gqu))
       gqu <- x$qu[which]
   if (missing(gth))
       gth <- quantile(x$gumbel[, which], gqu)
   dependent <- (1:(dim(x$data)[[2]]))[-which]
   if (length(gqu) < length(dependent))
       gqu <- rep(gqu, length = length(dependent))
   qfun <- function(X, yex, wh) {
       Q <- function(yex, ydep, param) {
           a <- param[1]
           b <- param[2]
           m <- param[3]
           s <- param[4]
           obj <- function(yex, ydep, a, b, m, s) {
               mu <- a * yex + m * yex^b
               sig <- s * yex^b
               log(sig) + 0.5 * ((ydep - mu)/sig)^2
           }
           res <- sum(obj(yex, ydep, a, b, m, s))
               if (is.infinite(res)){
                       if (res < 0){ res <- -(10^8) }
                       else res <- 10^8
                       warning("Infinite value of Q in cmvxDependence")
               }
           res
       }



       o <- try(optim(c(0.5, 0, 0, 1), Q, lower = c(10^(-8), -(10^8),
           -(10^8), 10^(-8)), upper = c(1, 10^8, 10^8, 10^8), method = "L-BFGS-B",
           yex = yex[wh], ydep = X[wh]))



       if (class(o) == "try-error" || o$convergence != 0) {
           warning("Non-convergence in cmvxDependence")
           o$par <- rep(NA, 4)
       }
       if (!is.na(o$par[1]))
           if (o$par[1] < 10^(-5) & o$par[2] < 0) {
               Q <- function(yex, ydep, param) {
                 param <- param[-1]
                 b <- param[1]
                 cee <- param[2]
                 d <- param[3]
                 m <- param[4]
                 s <- param[5]
                 obj <- function(yex, ydep, a, b, cee, d, m,
                   s) {
                   mu <- cee - d * log(yex) + m * yex^b
                   sig <- s * yex^b
                   log(sig) + 0.5 * ((ydep - mu)/sig)^2
                 }
                 res <- sum(obj(yex, ydep, a, b, cee, d, m,
                   s))
                 res
               }
               o <- try(optim(c(0, 0, 0, 0, 0, 1), Q, lower = c(10^(-8),
                 -Inf, -Inf, 10^(-8), -Inf, 10^(-8)), upper = c(1,
                 10^8, 10^8, 1 - 10^(-8), Inf, Inf), method = "L-BFGS-B",
                 yex = yex[wh], ydep = X[wh]))
               if (class(o) == "try-error" || o$convergence != 0) {
                 warning("Non-convergence in cmvxDependence")
                 o$par <- rep(NA, 4)
               }
           }
           else o$par <- c(o$par[1:2], 0, 0)
       o$par[1:4]
   }
   yex <- c(x$gumbel[, which])
   wh <- yex > gth
   res <- apply(x$gumbel[, dependent], 2, qfun, yex = yex, wh = wh)
   dimnames(res)[[1]] <- letters[1:4]
   gdata <- x$gumbel[wh, -which]
   tfun <- function(i, data, yex, a, b, cee, d) {
       data <- data[, i]
       a <- a[i]
       b <- b[i]
       cee <- cee[i]
       d <- d[i]
       if (is.na(a))
           rep(NA, length(data))
       else {
           if (a < 10^(-5) & b < 0)
               a <- cee - d * log(yex)
           else a <- a * yex
           (data - a)/(yex^b)
       }
   }
   z <- try(sapply(1:(dim(gdata)[[2]]), tfun, data = gdata,
       yex = yex[wh], a = res[1, ], b = res[2, ], cee = res[3,
           ], d = res[4, ]))
   if (class(z) %in% c("Error", "try-error")) {
       z <- matrix(nrow = 0, ncol = dim(x$data)[[2]] - 1)
   }
   else if (is.R()) {
       if (!is.array(z)) {
           z <- matrix(nrow = 0, ncol = dim(x$data)[[2]] - 1)
       }
   }
   else dimnames(z)[[2]] <- dimnames(res)[[2]]
   res <- list(call = theCall, parameters = res, Z = z, gth = gth,
       gqu = gqu, which = which)
   oldClass(res) <- "cmvxDependence"
   res
}

