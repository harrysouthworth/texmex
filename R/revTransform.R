`revTransform` <-
function (x, data, qu, th = 0, sigma = 1, xi = 0, method = "mixture") {
   if (!is.element(method, c("mixture", "empirical")))
       stop("method should be 'mixture' or 'empirical'")

   n <- length(data)
   probs <- (1:n)/(n + 1)
   px <- sapply(x, function(x, p) p[abs(x - p) == min(abs(x - p))][1], p = probs) # take 1st item in case of ties
   px <- as.integer(round(px * (1 + n)))
   res <- sort(data)[px]
   if (method == "mixture") {
     res[res > th] <- u2gpd(x[res > th], p=1-qu, th = th, sigma = sigma, xi = xi)
   }
   res[order(x)] <- sort(res)
   res
}

test.revTransform <- function(){
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
               
  checkEqualsNumeric(x.l.e,x,tol=0.0001,msg="revTransform: empirical transformation, laplace target")
  checkEqualsNumeric(x.l.m,x,tol=0.0001,msg="revTransform: mixture transformation, laplace target")
  checkEqualsNumeric(x.g.e,x,tol=0.0001,msg="revTransform: empirical transformation, gumbel target")
  checkEqualsNumeric(x.g.m,x,tol=0.0001,msg="revTransform: mixture transformation, gumbel target")
}

