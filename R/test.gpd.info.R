test.gpd.info <-
function(){
  lmod <- evm(ALT.M, data=liver, qu=.5, xi=~I(240*as.numeric(dose)),
              cov="numeric")
  checkTrue(all(sqrt(diag(solve(gpd.info(lmod)))) > 0),
            msg="gpd.inf: SDs positive")
  
  # Check equality to numerical approximation in big samples
  set.seed(20110923)
  tol <- 10^(-3)
  for (i in 1:10){
    x <- rt(10000, 10)
    junk <- evm(x, qu=.9, penalty="none", cov="numeric")
    msg <- paste("gpd.info: t", i, "equality to numerical", sep="")
    checkEqualsNumeric(junk$cov, solve(gpd.info(junk)), tolerance=tol,
                       msg=msg)
    
    # check estimation when we have a penalty
    gp1 <- list(c(0, 0), diag(c(10^4, .05)))
    gp2 <- list(c(0, 0), diag(c(.1, 10^4)))
    junk1 <- evm(x, qu=.9, priorParameters = gp1, cov="numeric")
    junk2 <- evm(x, qu=.9, priorParameters = gp2, cov="numeric")
    msg1 <- paste("gpd.info: t", i, "equality to numerical, penalty on xi", sep="")
    msg2 <- paste("gpd.info: t", i, "equality to numerical, penalty on phi", sep="")
    tol <- 0.01
    checkEqualsNumeric(junk1$cov, solve(gpd.info(junk1)), tolerance=tol, msg=msg1)
    checkEqualsNumeric(junk2$cov, solve(gpd.info(junk2)),
                       tolerance=tol, msg=msg2)
    
    # check estimation when we have covariates
    n <- 10000
    x <- 1/runif(n)
    data <- data.frame(x=x,y=rexp(n,exp(2 + x)))
    
    junk3 <- evm(y,data=data,phi =~ x,th=0)
    msg3 <- paste("gpd.info: t",i,"equality to numerical, covariates in phi",sep="")
    checkEqualsNumeric(junk3$cov, solve(gpd.info(junk3)), tolerance=tol, msg=msg3)
    
    x <- runif(n,-0.5,0.5)
    data <- data.frame(x=x,y = rgpd(n,sigma = exp(3+2*x), xi=x))
    
    junk4 <- evm(y,data=data,phi=~x, xi = ~ x,th=0)
    msg4 <- paste("gpd.info: t",i,"equality to numerical, covariates in phi and xi",
                  sep="")
    checkEqualsNumeric(junk4$cov, solve(gpd.info(junk4)), tolerance=tol, msg=msg4)
  }
}
