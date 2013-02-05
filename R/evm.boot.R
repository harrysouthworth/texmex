evmBoot <- function(o, R=1000, trace=100, theCall){
    if (class(o) != "evmOpt"){
        stop("o must be of class 'evmOpt'")
    }

    if (missing(theCall)){ theCall <- match.call() }

    D <- o$data$D
    param <- texmexGetParam(D, o$coefficients)

    param <- lapply(1:length(D), function(i, co, d){
                                           d[[i]] %*% co[[i]]
                                       }, co=param, d=D)
    param <- do.call("cbind", param)

    rng <- o$family$rng

    bfun <- function(i, param, o){
        if (i %% trace == 0){ cat("Replicate", i, "\n") }
        d <- o$data
        d$y <- rng(nrow(param), param, o)
        
		evm.fit(d, o$family, penalty=o$penalty,
		        priorParameters=o$priorParameters,
		        start=o$coefficients, hessian=FALSE)$par
    }

    res <- t(sapply(1:R, bfun, param=param, o=o))

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
	}
	invisible()
}

show.evmBoot <- print.evmBoot
show.summary.evmBoot <- print.summary.evmBoot

test.evmBoot <- function(){
    set.seed(20111007)
    # Compare bootstrap standard errors with those given by Coles
    # page 85
    tol <- 0.1
    cse <- c(.958432, .101151)
    raingpd <- gpd(rain, th=30, penalty="none")
    rainboot <- bootgpd(raingpd, R=100, trace=100)
    rainrep <- rainboot$replicates
    rainrep[,1] <- exp(rainrep[, 1])
    bse <- apply(rainrep, 2, sd)

    checkEqualsNumeric(cse[1],bse[1],tolerance=tol,msg="bootgpd: rain se(sigma) matches Coles")
    checkEqualsNumeric(cse[2],bse[2],tolerance=tol,msg="bootgpd: rain se(xi) matches Coles")

    # Check bootstrap medians are close to point estimates (the MLEs are
    # biased and the distribution of sigma in particular is skewed, so use
    # medians, not means, and allow a little leeway

    best <- apply(rainrep, 2, median)
    cest <- coef(raingpd); cest[1] <- exp(cest[1])
    checkEqualsNumeric(cest[1],best[1],tolerance=tol,msg="bootgpd: rain median of sigma matches point estimate")
    checkEqualsNumeric(cest[2],best[2],tolerance=tol,msg="bootgpd: rain medians of xi matches point estimate")

    ##################################################################
    # Do some checks for models with covariates. Due to apparent instability
    # of the Hessian in some cases, allow some leeway

    lmod <- gpd(log(ALT.M / ALT.B), data=liver, qu=.7,
                xi= ~ as.numeric(dose), phi= ~ as.numeric(dose))
    lboot <- bootgpd(lmod, R=200, trace=100)
    bse <- apply(lboot$replicates, 2, sd)
    rse <- bse / lmod$se
    rse <- ifelse(rse < 1, 1/rse, rse)
    checkTrue(max(rse) < 1.5, msg="bootgpd: SEs with xi in model")

    best <- apply(lboot$replicates, 2, median)
    rest <- best / coef(lmod)
    rest <- ifelse(rest < 1, 1/rest, rest)
    checkTrue(max(rest) < 1.5, msg="bootgpd: medians in line with point ests")

    ## Check penalization works - set harsh penalty and do similar
    ## checks to above

    pp <- list(c(0, .5), diag(c(.5, .05)))
    raingpd <- gpd(rain, th=30, penalty="none", priorParameters=pp)
    rainboot <- bootgpd(raingpd, R=1000, trace=100)

    bse <- apply(rainboot$replicates, 2, sd)
    rse <- bse / raingpd$se
    rse <- ifelse(rse < 1, 1/rse, rse)
    checkTrue(max(rse) < 1.1, msg="bootgpd: SEs with xi in model")

    best <- apply(rainboot$replicates, 2, median)
    rest <- best / coef(raingpd)
    rest <- ifelse(rest < 1, 1/rest, rest)
    checkTrue(max(rest) < 1.1, msg="bootgpd: medians in line with point ests")
}

