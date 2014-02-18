test.copula <-
function(){
  fun <- function(d) apply(d,2,function(x)(1:n)[rank(x)])/(1+n)
  n <- 200
  
  u2 <- cbind(sample(n),sample(n))
  d2 <- fun(u2)
  
  u3 <- cbind(sample(n),sample(n),sample(n))
  d3 <- fun(u3)
  
  checkEqualsNumeric(d2,copula(u2)$copula,msg="copula: 2 dimensional")
  checkEqualsNumeric(d3,copula(u3)$copula,msg="copula: 3 dimensional")
  
  op <- options()
  options(show.error.messages=FALSE)
  checkException(copula(TRUE),msg="copula: exception")
  checkException(copula("text"),msg="copula: exception")
  options(op)
}
