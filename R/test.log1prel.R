test.log1prel <-
function() {
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
