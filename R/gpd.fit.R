gpd.fit <- function(y, th, X.phi, X.xi, penalty="none", start=NULL,
        priorParameters = NULL, maxit = 10000, trace = 0) {

    gpd.lik <- function(par, y, th, X.phi, X.xi, penalty = "none", 
        priorParameters = NULL) {
        keepsc <- par[1:ncol(X.phi)]
        keepxi <- par[-(1:ncol(X.phi))]
        sc <- colSums(par[1:ncol(X.phi)] * t(X.phi))
        xi <- colSums(par[(ncol(X.phi) + 1):(ncol(X.phi) + ncol(X.xi))] * 
            t(X.xi))
        y <- (y - th)/exp(sc)
        y <- 1 + xi * y
        if (min(y) <= 0) {
            l <- 10^6
        }
        else {
            l <- sum(sc) + sum(log(y) * (1/xi + 1))
        }
        if (penalty == "none") {
            l <- l
        }
        else if (penalty %in% c("quadratic", "gaussian")) {
            p <- mahalanobis(matrix(c(keepsc, keepxi), nrow = 1), 
                center = priorParameters[[1]], cov = priorParameters[[2]])
            l <- l + p
        }
        else if (penalty %in% c("lasso", "l1", "laplace")) {
            p <- sum(abs(c(keepsc, keepxi) - priorParameters[[1]]) * 
                diag(priorParameters[[2]]))
            l <- l + p
        }
        else stop("penalty can be 'none', 'lasso' or 'gaussian'")
        l
    }

    if (is.null(start)){
        start <- c(log(mean(y)), rep(1e-05, -1 + ncol(X.phi) + ncol(X.xi)))
    }

    o <- optim(par = start, fn = gpd.lik, y = y, X.phi = X.phi, 
        X.xi = X.xi, th = th, penalty = penalty, control = list(maxit = maxit, 
            trace = trace), priorParameters = priorParameters, 
        hessian = TRUE)
    invisible(o)
}
