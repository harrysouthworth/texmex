##' Miscellaneous internal functions
NULL

##' Accurately compute (exp(x) - 1) / x
##' @param x numeric vector
##' @return numeric vector
.exprel <- function(x) {
  .Call(.c.exprel, x)
}

##' Accurately compute log(1 + x) / x
##' @param x numeric vector
##' @return numeric vector
.log1prel <- function(x) {
  .Call(.c.log1prel, x)
}

##' Compute pmax(x y, -1) in such a way that zeros in x beat
##' infinities in y.
##'
##' This is a common pattern in much of the distribution code, so it's
##' worth factoring out.
##' @param x a numeric vector
##' @param y a numeric vector
##' @return an appropriate numeric vector
.specfun.safe.product <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  xy <- x * y
  xy[(x==0) & is.infinite(y)] <- 0
  pmax(xy, -1)
}


test.exprel <- function(x) {
  ## pull in from the namespace because RUnit doesn't give access to
  ## internal functions.
  exprel <- texmex:::.exprel

  ## first check some simple values
  values <- c(-Inf, NA, seq(-5, 5, length.out=10))
  r.values <- expm1(values) / values
  exprel.values <- exprel(values)
  checkEqualsNumeric(exprel.values, r.values,
                     msg="exprel: value tests")


  ## the naive R code doesn't compute these right so force the
  ## computation.
  values <- c(Inf, 0)
  exprel.values <- exprel(values)
  desired.values <- c(Inf, 1)
  checkEqualsNumeric(exprel.values, desired.values,
                     msg="exprel: special value tests")

}


test.log1prel <- function() {
  ## pull in from the namespace

  log1prel <- texmex:::.log1prel
  exprel   <- texmex:::.exprel

  ## check simple values

  x <- runif(10, -1, 5)
  log1prel.values <- log1prel(x)

  r.values <- rep.int(1, length(x))
  flag <- x!=0
  r.values[flag] <- log1p(x[flag]) / x[flag]
  checkEqualsNumeric(log1prel.values, r.values,
                     msg="log1prel: value tests")

  ## special values
  checkEqualsNumeric(log1prel(c(-1, 0, Inf)), c(Inf, 1, 0),
                     msg="log1prel: special value tests")

  ## check with reference to .exprel

  x <- runif(10, -5, 5)
  y <- runif(length(x), -5, 5)

  xi <- exprel(x * y) * y

  checkEqualsNumeric(log1prel(xi * x) * xi,
                     y, "log1prel: exprel inversion")

}


test.specfun.safe.product <- function() {
  prod <- texmex:::.specfun.safe.product

  ## simple values
  x <- runif(10, -5, 5)
  y <- runif(length(x), -5, 5)

  checkEqualsNumeric(prod(x, y), pmax(x*y, -1),
                     "safe product: simple values")

  ## complicated values
  x <- c(0, 0, 1, 1)
  y <- c(Inf, -Inf, Inf, -Inf)

  res <- c(0, 0, Inf, -1)
  checkEqualsNumeric(prod(x, y), res,
                     "safe product: complicated values")
}
