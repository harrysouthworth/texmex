test.MCS <-
function(){
  myMCS <- function(x,p){
    # First and second args are
    # x (dxn matrix) and p (vector of probabilities).
    
    n <- dim(x)[2]
    d <- dim(x)[1]
    u <- t(apply(x,1,rank))/(n+1)# schmid and schmidt use n not n+1
    
    rho <- numeric(length(p))
    for(k in 1:length(p)) {
      Diff <- p[k] - u
      Diff[Diff<0] <- 0
      Prod <- apply(Diff,2,prod)
      
      num <- mean(Prod) - ((p[k]^2)/2)^d
      den <- (p[k]^(d+1)) / (d+1) - (0.5*p[k]^2)^d
      rho[k] <- num/den
    }
    
    rho
  }
  # simulated data - dimension 2
  n <- 1000
  by <- 0.01
  p <- seq(by,1-by,by=by)
  data <- rbind(rnorm(n),rnorm(n))
  
  tmRl <- MCS(t(data),p)
  myRl <- myMCS(data,tmRl$p)
  
  checkEqualsNumeric(myRl,tmRl$mcs,msg="MCS: independent normal data")
  checkEqualsNumeric(p,tmRl$p, msg="MCS: mathching p argument")
  
  # winter air pollution data - dimension 5
  tmWinterMCS <- MCS(winter,p)
  myWinterMCS <- myMCS(t(winter),p)
  checkEqualsNumeric(myWinterMCS,tmWinterMCS$mcs,msg="MCS: winter air pollution data")
  
  # summer airpollution data - dimension 5
  tmSummerMCS <- MCS(summer,p)
  mySummerMCS <- myMCS(t(summer),p)
  checkEqualsNumeric(mySummerMCS,tmSummerMCS$mcs,msg="MCS: summer air pollution data")
}
