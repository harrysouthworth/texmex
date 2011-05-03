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

        y <- rgpd(length(xi), xi=xi, sigma=exp(phi))
        gpd.fit(y=y,th=min(y) -1, X.phi=X.phi, X.xi=X.xi, penalty=prior, priorParameters=pp, start=co)$par
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

summary.bootgpd <- function(x, ...){
    print(x)
    cat("\nCorrelation:\n")
    print(cor(x$replicates))
    invisible(list(mean=apply(x$replicates, 2, mean), cov=cov(x$replicates)))
}



show.bootgpd <- print.bootgpd
