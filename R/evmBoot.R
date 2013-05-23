evmBoot <- function(o, R=1000, trace=100, theCall){
    if (class(o) != "evmOpt"){
        stop("o must be of class 'evmOpt'")
    }

    if (missing(theCall)){ theCall <- match.call() }

    d <- o$data
    param <- texmexMakeParams(coef(o), d$D)
    rng <- o$family$rng

    bfun <- function(i){
        if (i %% trace == 0){ cat("Replicate", i, "\n") }

        d$y <- rng(length(d$y), param, o)

        evmFit(d, o$family, th=o$threshold, prior=o$penalty,
               priorParameters=o$priorParameters,
               start=o$coefficients, hessian=FALSE)$par
    }

    res <- t(sapply(1:R, bfun))

    se <- apply(res, 2, sd)
    b <- apply(res, 2, mean) - coef(o)

    if (any(abs(b/se) > .25)){
        warning("Ratio of bias to standard error is high")
    }

    res <- list(call=theCall, replicates=res, map=o)
    oldClass(res) <- "evmBoot"
    res
}

print.evmBoot <- function(x, ...){
    print(x$call)
    means <- apply(x$replicates, 2, mean)
    medians <- apply(x$replicates, 2, median)
    sds <- apply(x$replicates, 2, sd)
    bias <- means - x$map$coefficients
    res <- rbind(x$map$coefficients, means, bias, sds, medians)
    rownames(res) <- c("Original", "Bootstrap mean", "Bias", "SD", "Bootstrap median")
    #colnames(res) <- names(summary(rnorm(3)))
    print(res, ...)
    if (any(abs(res[3,] / res[4,]) > .25)){
        warning("Ratio of bias to standard error is high")
    }
    invisible(res)
}


coefficients.evmBoot <- coef.evmBoot <- function(object, ...){
    apply(object$replicates, 2, mean)
}

summary.evmBoot <- function(object, ...){
    means <- apply(object$replicates, 2, mean)
    medians <- apply(object$replicates, 2, median)
    sds <- apply(object$replicates, 2, sd)
    bias <- means - coef(object$map)
    res <- rbind(coef(object$map), means, bias, sds, medians)
    rownames(res) <- c("Original", "Bootstrap mean", "Bias", "SD", "Bootstrap median")

    if (any(abs(res[3,] / res[4,]) > .25)){
        warning("Ratio of bias to standard error is high")
    }

	covs <- var(object$replicates)
	res <- list(call = object$call, margins=res, covariance=covs)
	oldClass(res) <- "summary.evm.boot"
    res
}

print.summary.evmBoot <- function(x, ...){
	print(x$call)
	print(x$margins)
	cat("\nCorrelation:\n")
    print(cov2cor(x$covariance))
    invisible()
}

plot.evmBoot <- function(x, col=4, border=FALSE, ...){
	pfun <- function(x, col, border, xlab,...){
		d <- density(x, n=100)
		hist(x, prob=TRUE, col=col, border=border, main="", xlab=xlab, ...)
		lines(d, lwd=2, col="grey")
		rug(x)
		invisible()
	}
	for (i in 1:ncol(x$replicates)){
		pfun(x$replicates[,i], xlab=colnames(x$replicates)[i], col=col, border=border)
		abline(v=coef(x$map)[i], col="cyan")
	}
	invisible()
}

show.evmBoot <- print.evmBoot
show.summary.evmBoot <- print.summary.evmBoot

test.evmBoot <- function(){
  tol <- 0.1
  
  for(Family in list(gpd,gev)){
    set.seed(20111007)

    pst <- function(msg) texmexPst(msg,Family=Family)

    u    <- switch(Family$name,GPD=30,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    fit <- evm(data,th=u,family=Family, penalty="none")
    boot <- evmBoot(fit, R=200, trace=1000)
    co <- coef(fit)
    rep <- boot$replicates
    scaleColumn <- switch(Family$name,GPD=1,GEV=2)
    rep[,scaleColumn] <- exp(rep[, scaleColumn])
    
    # Compare bootstrap standard errors with those given by Coles
    # pages 59 and 85 for GEV and GPD resp
    
    bse <- apply(rep, 2, sd)
    cse <- switch(Family$name,GPD=c(.958432, .101151),GEV=c(0.02792848, 0.02024846, 0.09823441))
    
    checkEqualsNumeric(cse,bse,tolerance=tol,
                       msg=pst("evmBoot: bootstrap se of parameter estimates matches Coles"))

    ## Check penalization works - set harsh penalty and do similar
    ## checks to above

    pp <- switch(Family$name,GPD=list(c(0, .5), diag(c(.5, .05))), GEV=list(c(5, 0, .5), diag(c(.5, .5, .05))))
    fit <- evm(data, th=u, penalty="none", priorParameters=pp,family=Family)
    boot <- evmBoot(fit, R=1000, trace=1100)
    
    bse <- apply(boot$replicates, 2, sd)
    rse <- bse / fit$se
    rse <- ifelse(rse < 1, 1/rse, rse)
    checkTrue(max(rse) < 1.1, msg=pst("evmBoot: SEs with xi in model, with penalty applied"))
    
    best <- apply(boot$replicates, 2, median)
    rest <- best / coef(fit)
    rest <- ifelse(rest < 1, 1/rest, rest)
    checkTrue(all(rest < switch(Family$name,GPD=c(1.1,1.1),GEV=c(1.1,1.1,1.3))), msg=pst("evmBoot: medians in line with point ests, with penalty applied"))
    
    ##################################################################
    # models with covariates. Due to apparent instability
    # of the Hessian in some cases, allow some leeway

    n <- 1000
    mu <- 1
    phi <- 5
    xi <- 0.05
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    
    test <- function(boot,fit,txt){
      bse <- apply(boot$replicates, 2, sd)
      rse <- bse / fit$se
      rse <- ifelse(rse < 1, 1/rse, rse)
      checkTrue(max(rse) < 1.5, msg=pst(paste("evmBoot: SEs with covariates in",txt)))
    
      best <- apply(boot$replicates, 2, median)
      rest <- best / coef(fit)
      rest <- ifelse(rest < 1, 1/rest, rest)
      checkTrue(all(rest < switch(Family$name,GPD=1.5,GEV=c(1.5,1.5,4,1.5))), msg=pst(paste("evmBoot: medians in line with point ests, covariates in",txt)))
    }
    
    param <- switch(Family$name,GPD=cbind(X[,1],xi),GEV=cbind(mu,X[,1],xi))
    start <- switch(Family$name,GPD=c(0,1,0.05),GEV=c(1,0,1,0.05))
    X$Y <- Family$rng(n,param,list(threshold=th))
    
    fit <- evm(Y,data=X,phi=~a,th=th,family=Family,start=start)
    boot <- evmBoot(fit, R=200, trace=201)
  
    test(boot,fit,"phi")
    
    param <- switch(Family$name,GPD=cbind(phi,X[,2]),GEV=cbind(mu,phi,X[,2]))

    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,xi=~b,th=th,family=Family)
    boot <- evmBoot(fit, R=200, trace=201)
    test(boot,fit,"xi")
  }
}


