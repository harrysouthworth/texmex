bootgpd <- function(x, R=99, trace=10){
    if (class(x) != "gpd"){
        stop("x must be of class 'gpd'")
    }
    
    theCall <- match.call()

    nphi <- ncol(x$X.phi); nxi <- ncol(x$X.xi)
    phi <- coef(x)[1:nphi]
    xi <- coef(x)[(nphi+1):(nphi+nxi)]
    phi <- colSums(phi * t(x$X.phi))
    xi <- colSums(xi *t(x$X.xi))

    if ("priorParameters" %in% names(x$call)){
        pp <- x$priorParameters
    }
    else {
        pp <- list(rep(0, nphi+nxi), diag(rep(10^4, nphi+nxi)))
    }

    bfun <- function(i, xi, phi, X.phi, X.xi, co, pp, prior){
        if (i %% trace == 0){ cat("Replicate", i, "\n") }

        r <- rgpd(length(xi), xi=xi, sigma=exp(phi))
		wh <- gpd.fit(y=r,th=min(r), X.phi=X.phi, X.xi=X.xi, penalty=prior, priorParameters=pp, start=co)$par
    }

    res <- t(sapply(1:R, bfun, xi=xi, phi=phi, X.phi=x$X.phi, X.xi=x$X.xi, co=coef(x), pp=pp, prior=x$penalty))

    res <- list(call=theCall, replicates=res)

    oldClass(res) <- "bootgpd"
    res
}

print.bootgpd <- function(x, ...){
    print(x$call)
    res <- t( apply(x$replicates, 2, summary))
    #colnames(res) <- names(summary(rnorm(3)))
    print(res, ...)
    invisible(res)
}

summary.bootgpd <- function(object, ...){
    mars <- t(apply(object$replicates,2, summary))
	covs <- var(object$replicates)
	res <- list(call = object$call, margins=mars, covariance=covs)
	oldClass(res) <- "summary.bootgpd"
    res
}

print.summary.bootgpd <- function(x, ...){
	print(x$call)
	print(x$margins)
	cat("\nCorrelation:\n")
    print(cov2cor(x$covariance))
    invisible()
}

plot.bootgpd <- function(x, col="blue", border=FALSE, ...){
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

show.bootgpd <- print.bootgpd
show.summary.bootgpd <- print.summary.bootgpd