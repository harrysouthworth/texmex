gpd.fit <- function(y, th, X.phi, X.xi, penalty="none", start=NULL,
        priorParameters = NULL, maxit = 10000, trace = 0, scale.=TRUE) {

    if (scale.){
        # Scale and centre the data in the hope of stabilizing computations
        whphi <- (1:ncol(X.phi))[colnames(X.phi) != "(Intercept)"]
        whxi <- (1:ncol(X.xi))[colnames(X.xi) != "(Intercept)"]
        sphi <-apply(X.phi, 2, sd)
        sxi <- apply(X.xi, 2, sd)
        mphi <- apply(X.phi, 2, mean)
        mxi <- apply(X.xi, 2, mean)
        X.phi[,whphi] <- apply(matrix(X.phi[,whphi]), 2, scale)
        X.xi[,whxi] <- apply(matrix(X.xi[,whxi]), 2, scale)
    }
    
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

    if (scale.){
        # Need to account for scaling of the data...
        # First do the parameter estimates...

        phiPar <- o$par[1:ncol(X.phi)]
        xiPar <- o$par[(ncol(X.phi)+1):(ncol(X.phi) + ncol(X.xi))]

        phiPar[whphi] <- phiPar[whphi] / sphi[c(1:ncol(X.phi))[whphi]]
        if (colnames(X.phi)[1] == "(Intercept)"){
            phiPar[1] <- phiPar[1] - sum(phiPar[-1] * mphi[(1:ncol(X.phi))[whphi]])
        }
        xiPar[whxi] <- xiPar[whxi] / sxi[(1:ncol(X.xi))[whxi]]
        if (colnames(X.xi)[1] == "(Intercept)"){
            xiPar[1] <- xiPar[1] - sum(xiPar[-1] * mxi[(1:ncol(X.xi))[whxi]])
        }

        o$par <- c(phiPar, xiPar)

        ### Scale the covariance...
        U <- solve(o$hessian)

        uphi <- matrix(U[1:ncol(X.phi), 1:ncol(X.phi)], ncol=ncol(X.phi))
        ssphi <-  c(1 + (uphi[1] + sum(colSums(uphi[-1,-1])) - 2*sum(uphi[1,])) / uphi[1,1], 1/sphi[-1])
browser()

        uxi <- matrix(U[(ncol(X.phi)+1):(ncol(X.phi) + ncol(X.xi)), (ncol(X.phi)+1):(ncol(X.phi) + ncol(X.xi))], ncol=ncol(X.xi))
        ssxi <-  c(1 + (uxi[1] + sum(colSums(uxi[-1,-1])) - 2*sum(uxi[1,])) / uxi[1,1], 1/sxi[-1])
        Dstar <- diag(c(ssphi, ssxi))
        S <- Dstar %*% U %*% Dstar
        o$hessian <- solve(S)

    } # Close if (scale.


    invisible(o)
}
