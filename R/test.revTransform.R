context("revTransform")

test_that("revTransform behaves as it should", {
    set.seed(20111010)
  n <- 5000
  x <- cbind(rexp(n),rexp(n,3))
  
  x.fit <- migpd(x,mqu = 0.5,penalty="none")
  
  y.l.m <- mexTransform(x.fit,method="mixture", margins="laplace")$transformed
  y.l.e <- mexTransform(x.fit,method="empirical",margins="laplace")$transformed
  
  y.g.m <- mexTransform(x.fit,method="mixture", margins="gumbel")$transformed
  y.g.e <- mexTransform(x.fit,method="empirical",margins="gumbel")$transformed
  
  distFun.l <- function(x) ifelse(x<0, exp(x)/2, 1-exp(-x)/2)
  distFun.g <- function(x) exp(-exp(-x))
  
  u.g.m <- distFun.g(y.g.m)
  u.g.e <- distFun.g(y.g.e)
  
  u.l.m <- distFun.l(y.l.m)
  u.l.e <- distFun.l(y.l.e)
  
  x.l.m <- cbind(revTransform(u.l.m[,1],x[,1],x.fit$mqu[1],x.fit$mth[1],exp(x.fit$models[[1]]$coefficients[1]), x.fit$models[[1]]$coefficients[2]),
                 revTransform(u.l.m[,2],x[,2],x.fit$mqu[2],x.fit$mth[2],exp(x.fit$models[[2]]$coefficients[1]), x.fit$models[[2]]$coefficients[2]))
  x.l.e <- cbind(revTransform(u.l.e[,1],x[,1],method="empirical"),
                 revTransform(u.l.e[,2],x[,2],method="empirical"))
  
  x.g.m <- cbind(revTransform(u.g.m[,1],x[,1],x.fit$mqu[1],x.fit$mth[1],exp(x.fit$models[[1]]$coefficients[1]), x.fit$models[[1]]$coefficients[2]),
                 revTransform(u.g.m[,2],x[,2],x.fit$mqu[2],x.fit$mth[2],exp(x.fit$models[[2]]$coefficients[1]), x.fit$models[[2]]$coefficients[2]))
  x.g.e <- cbind(revTransform(u.g.e[,1],x[,1],method="empirical"),
                 revTransform(u.g.e[,2],x[,2],method="empirical"))
  
  expect_that(x.l.e, equals(x),   expect_that(x.l.m, equals(x),   expect_that(x.g.e, equals(x),   expect_that(x.g.m, equals(x), }
)
