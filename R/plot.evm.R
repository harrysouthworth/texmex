plot.evmOpt <-
function(x, main=rep(NULL,4), xlab=rep(NULL,4), nsim=1000, alpha=.05, ...){
    if (!missing(main)){
        if (length(main) != 1 & length(main) != 4){
            stop("main should have length 1 or 4")
        }
        else if (length(main) == 1){ main <- rep(main, 4) }
    }

    if (all(sapply(x$data$D,ncol) == 1)){
        plot(ppevm(x, nsim=nsim, alpha=alpha), main=main[1], xlab=xlab[1])
        plot(qqevm(x, nsim=nsim, alpha=alpha), main=main[2], xlab=xlab[2])
        plotrl.evmOpt(x, main=main[3], xlab=xlab[3], smooth=FALSE, ...)
        plot(hist.evmOpt(x, main=main[4], xlab=xlab[4]))
    } else { # Covariates in the model
    
        np <- length(x$data$D)
        lp <- predict(x,type="lp", unique.=FALSE)
        Which <- as.logical(apply(lp[,1:np],2,var)) # identifies which cols have covariates

        x$data$y <- resid(x)
        x$threshold <- 0
        x$coefficients <- rep(0, length(x$data$D)) # phi not sigma, so 0 not 1
        plot(ppevm(x, nsim=nsim, alpha=alpha), main=main[1], xlab=xlab[1])
        plot(qqevm(x, nsim=nsim, alpha=alpha), main=main[2], xlab=xlab[2])

        
        for(i in (1:length(x$data$D))[Which]){
          ParName <- names(x$data$D[i])
          plot(lp[,i],resid(x),main=paste("Residuals vs fitted",ParName),xlab=paste("Fitted",ParName),ylab="Residuals")
          panel.smooth(lp[,i], resid(x), col.smooth=2)
        }
    }
    
    invisible()
}

test.plot.evmOpt <- function(){
  par(mfrow=c(2,2))
  mod <- evm(rain, th=30, penalty="none")
  res <- plot(mod,main=paste(rep("Figure 4.5 of Coles (2001)",4),
              c("\nProbability plot","\nQuantile Plot","\nReturn Level Plot\n(SCALE IS DAYS NOT YEARS)","\nDensity Plot")), 
              RetPeriodRange=c(3.65,365*10000))
  checkEquals(res,NULL,msg="plot.evmOpt: GPD successful execution")

  # check for very short tailed data
  set.seed(6)
  temp <- rgpd(1000,sigma=1,xi=-0.45)
  fit <- evm(temp,th=0)
  res <- plot(fit,main=c("GPD: PP","GPD: QQ","GPD: RL","GPD: Hist, Short tailed data"))
  checkEquals(res,NULL,msg="plot.evmOpt: GPD successful execution, short tailed data")
  
  # check for covariates in the model
  # GPD
  n <- 1000
  sig <- 2
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,sig,X[,2])
  X$Y <- Y
  fit <- evm(Y,data=X,xi=~b,th=0)
  res <- plot(fit)
  checkEquals(res,NULL,msg="plot.evmOpt: GPD with covariates successful execution")
 
  #GEV 
  # no covariates
  n <- 1000
  Y <- rgev(n,1,1,-.1)
  fit <- evm(Y,family=gev)
  par(mfrow=c(2,2))
  res <- plot(fit,main="GEV no covariates, neg xi")
  checkEquals(res,NULL,msg="plot.evmOpt: GEV no covariates, neg xi successful execution")
    
  Y <- rgev(n,1,1,.2)
  fit <- evm(Y,family=gev)
  res <- plot(fit,main="GEV no covariates, pos xi")
  checkEquals(res,NULL,msg="plot.evmOpt: GEV no covariates, pos xi successful execution")

  #GEV with covariates
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3),C= runif(n))
  Y <- rgev(n,X[,3],sig,X[,2])
  X$Y <- Y
  fit <- evm(Y,data=X,xi=~b,mu=~C,family=gev)
  par(mfrow=c(2,2))
  res <- plot(fit,main=rep("GEV with covariates",4))
  checkEquals(res,NULL,msg="plot.evmOpt: GEV with covariates successful execution")

  fit <- evm(Y,data=X,xi=~b,family=gev)
  par(mfrow=c(2,2))
  res <- plot(fit,main=rep("GEV with one covariate",4))
  checkEquals(res,NULL,msg="plot.evmOpt: GEV with one covariate successful execution")
}
