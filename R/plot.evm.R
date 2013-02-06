plot.evmOpt <-
function(x, main=rep(NULL,4), xlab=rep(NULL,4), nsim=1000, alpha=.05, ...){
    if (!missing(main)){
        if (length(main) != 1 & length(main) != 4){
            stop("main should have length 1 or 4")
        }
        else if (length(main) == 1){ main <- rep(main, 4) }
    }

    if (ncol(x$data$D$phi) == 1 && ncol(x$data$D$xi) == 1){
        ppevm(x, main=main[1], xlab=xlab[1], nsim=nsim, alpha=alpha)
        qqevm(x, main=main[2], xlab=xlab[2], nsim=nsim, alpha=alpha)
        plotrl.evmOpt(x, main=main[3], xlab=xlab[3], ...)
        hist.evmOpt(x, main=main[4], xlab=xlab[4])
    }
    else { # Covariates in the model
        fittedScale <- exp(coef(x)[1:ncol(x$data$D$phi)] %*% t(x$data$D$phi))
        fittedShape <- coef(x)[(ncol(x$data$D$phi) + 1):length(coef(x))] %*% t(x$data$D$xi)

        x$data$y <- resid(x)
        x$threshold <- 0
        x$coefficients <- rep(0, length(x$data$D)) # phi not sigma, so 0 not 1
        ppevm(x, main=main[1], xlab=xlab[1], nsim=nsim, alpha=alpha)
        qqevm(x, main=main[2], xlab=xlab[2], nsim=nsim, alpha=alpha)

        if(ncol(x$data$D[[1]]) > 1){ # phi
          plot(fittedScale,resid(x),
               main="Residuals vs Fitted Scale",
               xlab="Fitted scale",ylab="Residuals")
          panel.smooth(fittedScale, resid(x), col.smooth=2)
        }
        if(ncol(x$data$D[[2]]) > 1){ # xi
          plot(fittedShape,resid(x),
               main="Residuals vs Fitted Shape",
               xlab="Fitted shape",ylab="Residuals")
          panel.smooth(fittedShape, resid(x), col.smooth=2)
        }
    }

    invisible()
}

test.plot.evmOpt <- function(){
  par(mfrow=c(2,2))
  mod <- gpd(rain, th=30, penalty="none")
  res <- plot(mod,main=paste(rep("Figure 4.5 of Coles (2001)",4),
              c("\nProbability plot","\nQuantile Plot","\nReturn Level Plot\n(SCALE IS DAYS NOT YEARS)","\nDensity Plot")), 
              RetPeriodRange=c(3.65,365*10000))
  checkEquals(res,NULL,msg="plot.gpd: successful execution")
  
# check for very short tailed data
  set.seed(6)
  temp <- rgpd(1000,sigma=1,xi=-0.45)
  fit <- gpd(temp,th=0)
  res <- plot(fit,main=c("PP","QQ","RL","Hist, Short tailed data"))
  checkEquals(res,NULL,msg="plot.gpd: successful execution, short tailed data")
} 
