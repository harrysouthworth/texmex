# Author: Harry Southworth
# Date: 2011-11-25
## Purpose: Create a predict method for objects of class gpd and bgpd that
##          returns parameters, return levels or (maybe) return periods,
##          depending on arguments given.
#
# predict.gpd
# predict.bgpd
# predict.bootgpd
# rl
# rl.gpd
# rl.bgpd
# rl.bootgpd
# linearPredictors
# linearPredictors.gpd
# linearPredictors.bgpd
# linearPredictors.bootgpd

################################################################################
## gpd

predict.gpd <-
    # Get predictions for a gpd object. These can either be the linear predictors
    # or return levels.
function(object, M=1000, newdata=NULL, type="return level", se.fit=FALSE,
         ci.fit=FALSE, alpha=.050, unique.=TRUE, ...){
    theCall <- match.call()

    res <- switch(type,
                  "rl"=, "return level" = rl.gpd(object, M, newdata,
                                                 se.fit=se.fit, ci.fit=ci.fit,
                                                 alpha=alpha, unique.=unique.),
                  "lp" =,"link" = linearPredictors.gpd(object, newdata, se.fit,
                                                   ci.fit, alpha, unique.=unique.)
                  )
    res
}

## Linear predictor functions for GPD

linearPredictors.gpd <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                             alpha=.050, unique.=TRUE, full.cov=FALSE, ...){

    if (!is.null(newdata)){
        xi.fo <- object$formulae$xi
        phi.fo <- object$formulae$phi
        X.xi <-  if (!is.null(xi.fo)) { model.matrix(xi.fo,  newdata) } else { matrix(1, nrow(newdata)) }
        X.phi <- if (!is.null(phi.fo)){ model.matrix(phi.fo, newdata) } else { matrix(1, nrow(newdata)) }
    } else {
        X.xi <- object$X.xi
        X.phi <- object$X.phi
    }

    if (unique.){
        u <- !duplicated(cbind(X.phi,X.xi))
        X.xi  <- if(is.matrix(X.xi[u,]))  X.xi[u, ]  else if(ncol(X.xi) == 1)  cbind(X.xi[u,])  else t(cbind(X.xi[u,]))
        X.phi <- if(is.matrix(X.phi[u,])) X.phi[u, ] else if(ncol(X.phi) == 1) cbind(X.phi[u,]) else t(cbind(X.phi[u,]))
    }

    whichPhi <- 1:ncol(X.phi)
    whichXi  <- (ncol(X.phi) + 1):length(object$coefficients)
    phi <- c(object$coefficients[whichPhi] %*% t(X.phi))
    xi  <- c(object$coefficients[whichXi]  %*% t(X.xi))

    res <- cbind(phi, xi)

    if(ci.fit | se.fit | full.cov){
       phi.cov <- as.matrix(object$cov[whichPhi, whichPhi])
       xi.cov  <- as.matrix(object$cov[whichXi,  whichXi])

       if(ci.fit | se.fit){
          phi.se <- sqrt(rowSums((X.phi %*% phi.cov) * X.phi))
          xi.se <-  sqrt(rowSums((X.xi  %*% xi.cov)  * X.xi))
       }
    }

    if (ci.fit){
        z <- qnorm(1 - alpha/2)

        phi.lo <- phi - phi.se*z
        phi.hi <- phi + phi.se*z
        xi.lo <- xi - xi.se*z
        xi.hi <- xi + xi.se*z

        res <- cbind(res, phi.lo, phi.hi, xi.lo, xi.hi)
    } # Close if(ci.fit

    if (se.fit){
        res <- cbind(res, phi.se, xi.se)
    } # Close if(se.fit

    res <- addCov(res,X.phi)
    res <- addCov(res,X.xi)

    if (full.cov){ # Covariance of (phi, xi) for each unique (phi, xi) pair
        covar <- rep(0, nrow(X.xi))
        for(k in 1:length(covar)){
            for (i in 1:ncol(X.phi)){
                for (j in 1:ncol(X.xi)){
                    covar[k] <- covar[k] + X.phi[k, i] * X.xi[k, j] * object$cov[i, ncol(X.phi) + j]
                } # Close for j
            } # Close for i
        } # Close for k
        phi.var <- rowSums((X.phi %*% phi.cov) * X.phi)
        xi.var <- rowSums((X.xi %*% xi.cov) * X.xi)

        res <- list(link=res, cov=list(phi.var=phi.var, xi.var=xi.var, covariances=covar))
    }
    oldClass(res) <- "lp.gpd"
    res
}

## Return level functions for GPD

## Reversing arguments M and newdata for anyone who wants to call these functions
## directly

## Will want to get return levels when using GEV rather than GPD, so make
## rl generic

rl <- function(object, M = 1000, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, ...){
    UseMethod("rl")
}

linearPredictors <- function(object, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, ...){
    UseMethod("linearPredictors")
}

gpdDelta <- function(A, K){
   # This is not exact if a prior (penalty) function is used, but
   # the CI is approximate anyway.

    out <- matrix(0, nrow=3, ncol=length(K))

    if (A[3] == 0){ # exponential case
        out[1,] <- exp(A[2]) / A[1]
        out[2,] <- exp(A[2]) * log(K * A[1])
    } else {
        out[1,] <- exp(A[2]) * K^A[3] * A[1]^(A[3] - 1)
        out[2,] <- exp(A[2]) / A[3] * ((K*A[1])^A[3] - 1)
        out[3,] <- -exp(A[2]) / (A[3]*A[3]) * ( (K * A[1] )^A[3] - 1 ) +
                   exp(A[2]) / A[3] * (K * A[1])^A[3] * log(K * A[1])
    }

   out
}

rl.gpd <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                   alpha=.050, unique.=TRUE, ...){
    co <- linearPredictors.gpd(object, newdata=newdata, unique.=unique., full.cov=TRUE)
    covs <- co[[2]] # list(phi.var=phi.var, xi.var=xi.var, covariances=covar)
    co <- co[[1]]
    X <- co[,-(1:2)]
    if(is.null(dim(X))){
      X <- matrix(X)
      dimnames(X) <- list(dimnames(co)[[1]],dimnames(co)[[2]][-(1:2)])
    }

    gpdrl <- function(u, theta, phi, xi, m){
        res <- u + exp(phi) / xi *((m * theta)^xi -1)
        cbind(RL=res)
    }

    res <- lapply(M, gpdrl,
                  u=object$threshold, theta=object$rate, phi=co[,1], xi=co[,2])

    getse <- function(o, co, M){
        dxm <- lapply(split(co, 1:nrow(co)), gpdDelta, K=M)

        V <- lapply(1:nrow(co),
                    function(i, x, rate, n){
                        # Construct covariance matrix
                        cov <- matrix(c(x[[1]][i], rep(x[[3]][i], 2), x[[2]][i]), ncol=2)
                        matrix(c(rate * (1 - rate) / n, 0, 0,
                               0, cov[1,],
                               0, cov[2,]), ncol=3)
                    }, rate = o$rate, n = length(o$y) / o$rate, x=covs)

        # Get (4.15) of Coles, page 82, adjusted for phi = log(sigma)
        se <- sapply(1:length(V),
                     function(i, dxm, V){
                        V <- V[[i]]; dxm <- c(dxm[[i]])
                        sqrt(mahalanobis(dxm, center=c(0, 0, 0), cov=V, inverted=TRUE))
                     }, dxm=dxm, V=V)
        se
    }

    co <- cbind(rep(object$rate, nrow(co)), co)

    if (ci.fit){ # need to update plotrl.gpd too once profile lik confidence intervals implemented here
        ci.fun <- function(i, object, co, M, res, alpha){
            wh <- res[[i]];
            se <- getse(object, co, M[i])
            lo <- wh - qnorm(1 - alpha/2)*se
            hi <- wh + qnorm(1 - alpha/2)*se
            wh <- cbind(wh, lo=lo, hi=hi)

            colnames(wh) <- c("RL", paste(100*alpha/2, "%", sep = ""),
                              paste(100*(1 - alpha/2), "%", sep = ""))
            wh
        } # ci.fun
        res <- lapply(1:length(M), ci.fun, object=object, co=co, M=M, res=res, alpha=alpha)
    } # Close if (ci.fit

    if (se.fit){
        se.fun <- function(i, object, co, M, res, alpha){
            wh <- res[[i]]
            se <- getse(object, co, M[i])
            wh <- cbind(RL=wh, se.fit=se)
            wh
        } # ci.fun
        res <- lapply(1:length(M), se.fun, object=object, co=co, M=M, res=res, alpha=alpha)
    }

    cov.fun <- function(i,res){
      wh <- res[[i]]
      wh <- addCov(wh,X)
      wh
    }
    res <- lapply(1:length(M), cov.fun,res=res)

    names(res) <- paste("M.", M, sep = "")
    oldClass(res) <- "rl.gpd"
    res
}

################################################################################
## bgpd

predict.bgpd <- function(object, M=1000, newdata=NULL, type="return level",
                         se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE,
                         all=FALSE, sumfun=NULL, ...){
    theCall <- match.call()

    res <- switch(type,
                  "rl" = , "return level" = rl.bgpd(object, M=M, newdata=newdata,
                                                    se.fit=se.fit, ci.fit=ci.fit,
                                                    alpha=alpha, unique.=unique., all=all,
                                                    sumfun=sumfun,...),
                  "lp" = , "link" = linearPredictors.bgpd(object, newdata=newdata,
                                                      se.fit=se.fit, ci.fit=ci.fit,
                                                      alpha=alpha, unique.=unique., all=all,
                                                      sumfun=sumfun,...)
                  )
    res
}

linearPredictors.bgpd <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                              alpha=.050, unique.=TRUE, all=FALSE, sumfun=NULL, ...){
    if (!is.null(newdata)){
        xi.fo <- object$map$formulae$xi
        phi.fo <- object$map$formulae$phi
        X.xi <-  if (!is.null(xi.fo))  model.matrix(as.formula(xi.fo),  newdata) else matrix(1, nrow(newdata))
        X.phi <- if (!is.null(phi.fo)) model.matrix(as.formula(phi.fo), newdata) else matrix(1, nrow(newdata))
    } else {
        X.xi <- object$X.xi
        X.phi <- object$X.phi
    }

    X.phi.xi <- cbind(X.phi, X.xi)
    ModelHasCovs <- ncol(X.phi.xi) > 2
    if(ModelHasCovs){
      covCols <- apply(X.phi.xi, 2, function(x) !all(x==1))
      Xnames <- colnames(X.phi.xi)
      if(sum(covCols) == 1){
        X.phi.xi <- cbind(X.phi.xi[, covCols]) # cbind forces returned object to have a dim attribute
      } else {
        X.phi.xi <- X.phi.xi[, covCols]
      }
      colnames(X.phi.xi) <- Xnames[covCols]
    }

    if (unique.){
        X.phi.xi <- unique(X.phi.xi)
        u <- !duplicated(cbind(X.phi,X.xi))
        X.xi  <- if(is.matrix(X.xi[u,]))  X.xi[u, ]
                 else if(ncol(X.xi) == 1)  cbind(X.xi[u,])
                 else t(cbind(X.xi[u,]))
        X.phi <- if(is.matrix(X.phi[u,])) X.phi[u, ]
                 else if(ncol(X.phi) == 1) cbind(X.phi[u,])
                 else t(cbind(X.phi[u,]))
    }

    phi <- cbind(object$param[, 1:ncol(X.phi)])
    xi <- cbind(object$param[, (ncol(X.phi) + 1):(ncol(X.phi) + ncol(X.xi))])

    # Get point estimates (means)
    res <- lapply(1:nrow(X.xi),
                  function(i, phi, xi, Xphi, Xxi){
                      phi <- rowSums(t(t(phi) * c(Xphi[i, ])))
                      xi <- rowSums(t(t(xi) * c(Xxi[i, ])))
                      cbind(phi=phi, xi=xi)
                    }, phi=phi, xi=xi, Xphi=X.phi, Xxi=X.xi)

    ############################################################################
    ## Hard part should be done now. Just need to summarize

    if (ci.fit){
        if (is.null(sumfun)){
            sumfun <- function(x){
                c(mean(x), quantile(x, prob=c(.50, alpha/2,  1 - alpha/2)) )
            }
            neednames <- TRUE
        } else {
            neednames <- FALSE
        }

        res <- t(sapply(res, function(x, fun){ apply(x, 2, sumfun) }, fun=sumfun))

        if (neednames){
            nms <- c("Mean","50%", paste(100*alpha/2, "%", sep = ""),
                    paste(100*(1-alpha/2), "%", sep = ""))

            colnames(res) <- c(paste("phi:", nms), paste("xi:", nms))
        }
    } else if (se.fit){ warning("se.fit not implemented - ignoring") }
    else if (all){ res <- res }
    else { # Just point estimates
        res <- t(sapply(res, function(x){ apply(x, 2, mean) }))
    }

    if(!all){
      if(ModelHasCovs){
        res <- addCov(res,X.phi)
        res <- addCov(res,X.xi)
      }
    } else {
        if (ModelHasCovs & nrow(X.phi.xi) != length(res)){
            stop("Number of unique combinations of covariates doesn't match the number of parameters")
        }
        for (i in 1:length(res)){
          if(ModelHasCovs){
            res[[i]] <- cbind(res[[i]], matrix(rep(X.phi.xi[i,], nrow(res[[i]])), nrow=nrow(res[[i]]), byrow=TRUE))
            colnames(res[[i]]) <- c("phi", "xi", colnames(X.phi.xi))
          } else {
            colnames(res[[i]]) <- c("phi","xi")
          }
        }
    }

    oldClass(res) <- "lp.bgpd"
    res
}

rl.bgpd <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE, all=FALSE, sumfun=NULL,...){

    co <- linearPredictors.bgpd(object, newdata=newdata, unique.=unique., all=TRUE, sumfun=NULL)
    Covs <- linearPredictors.bgpd(object, newdata=newdata, unique.=unique., sumfun=NULL)
    X <- Covs[,-(1:2)]
    if(is.null(dim(X))){
      X <- matrix(X)
      dimnames(X) <- list(dimnames(Covs)[[1]],dimnames(Covs)[[2]][-(1:2)])
    }

    bgpdrl <- function(o, u, theta, m){
        res <- u + exp(o[, 1]) / o[, 2] *((m * theta)^o[, 2] -1)
        cbind(RL=res)
    }

    # co is a list with one element for each unique item in
    # new data. Need to loop over vector M and the elements of co

    getrl <- function(m, co, u, theta, ci.fit, alpha, all){
        res <- sapply(co, bgpdrl, u=u, theta=theta, m=m)

        if (ci.fit){
            if (is.null(sumfun)){
                sumfun <- function(x){
                    c(mean(x), quantile(x, prob=c(.50, alpha/2,  1 - alpha/2)))
                }
                neednames <- TRUE
            } else {
               neednames <- FALSE
            }
            res <- t(apply(res, 2, sumfun))
            if (neednames){
                colnames(res) <- c("Mean", "50%", paste(100*alpha/2, "%", sep = ""),
                                   paste(100*(1-alpha/2), "%", sep = ""))
            }

        } # Close if (ci.fit
        else if (se.fit){
            warning("se.fit not implemented")
            res <- apply(res, 2, mean)
        }
        else if (!all){ res <- apply(res, 2, mean) }
        res
    }

    res <- lapply(M, getrl, co=co, u=object$threshold, theta=object$map$rate, ci.fit=ci.fit, alpha=alpha, all=all)

    if(!all){
      cov.fun <- function(i,res){
        wh <- res[[i]]
        wh <- addCov(wh,X)
        wh
      }
      res <- lapply(1:length(M), cov.fun,res=res)
    }

    names(res) <- paste("M.", M, sep = "")
    oldClass(res) <- "rl.bgpd"
    res
}

################################################################################
## bootgpd

predict.bootgpd <- function(object, M=1000, newdata=NULL, type="return level",
                            se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE,
                            all=FALSE, sumfun=NULL, ...){
    theCall <- match.call()

    res <- switch(type,
                  "rl" = , "return level" = rl.bootgpd(object, newdata=newdata, M=M,
                                                       se.fit=se.fit, ci.fit=ci.fit,
                                                       alpha=alpha, unique.=TRUE,
                                                       all=all, sumfun=sumfun,...),
                  "lp" = , "link" = linearPredictors.bootgpd(object, newdata=newdata,
                                                         se.fit=se.fit, ci.fit=ci.fit,
                                                         alpha=alpha, unique.=TRUE,
                                                         all=all, sumfun=sumfun,...)
                  )
    res

}

addCov <- function(res,X){ # used in linearPredictors.* to add covariates to columns reported in output
  if(!is.null(dim(X))){
    if(dim(X)[2] > 1){
       cov <- X[,colnames(X) != "(Intercept)"]
       res <- cbind(res,cov)
       if(is.vector(cov)) colnames(res)[dim(res)[2]] <- colnames(X)[colnames(X) != "(Intercept)"]
    } else {
      if( any(X != 1) ){
        res <- cbind(res,X)
      }
    }
  }
  res
 }

namesBoot2bgpd <- function(bootobject){
    names(bootobject) <- c("call", "param", "original", "map")
    bootobject$X.phi <- bootobject$map$X.phi
    bootobject$X.xi <- bootobject$map$X.xi
    bootobject$threshold <- bootobject$map$threshold
    bootobject
}

linearPredictors.bootgpd <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=.050,
                                 unique.=TRUE, all=FALSE, sumfun=NULL,...){
    # This should just be the same as for a bgpd object, but some
    # names and stuff are different.
  object <- namesBoot2bgpd(object)
  res <- linearPredictors.bgpd(object, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit, all=all, unique.=unique., alpha=alpha, sumfun=sumfun,...)
  oldClass(res) <- "lp.bootgpd"
  res
}

rl.bootgpd <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=0.050, unique.=TRUE, all=FALSE, sumfun=NULL,...){
    # This should just be the same as for a bgpd object, but some
    # names are different.
  object <- namesBoot2bgpd(object)
  res <- rl.bgpd(object, M=M, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit,alpha=alpha, unique.=unique., all=all, sumfun=sumfun,...)
  oldClass(res) <- "rl.bootgpd"
  res
}

################################################################################
## Method functions

print.rl.gpd <- function(x, digits=3, ...){
    nms <- names(x)
    newnms <- paste("M =", substring(nms, 3), "predicted return level:\n")
    lapply(1:length(x), function(i, o, title){
                                 cat(title[i])
                                 print(o[[i]], digits=digits,...)
                                 cat("\n")
                                 NULL}, o=x, title=newnms)
    invisible(x)
}

summary.rl.gpd <- function(object, digits=3, ...){
    print.rl.gpd(object, digits=digits, ...)
}

print.rl.bgpd    <- print.rl.gpd
print.rl.bootgpd <- print.rl.gpd

summary.rl.bgpd    <- summary.rl.gpd
summary.rl.bootgpd <- summary.rl.gpd


print.lp.gpd <- function(x, digits=3, ...){
    cat("Linear predictors:\n")
    print(unclass(x), digits=3,...)
    invisible(x)
}

summary.lp.gpd <- function(object, digits=3, ...){
    print.lp.gpd(object, digits=3, ...)
}

summary.lp.gpd

summary.lp.bgpd    <- summary.lp.gpd
summary.lp.bootgpd <- summary.lp.gpd

print.lp.bgpd    <- print.lp.gpd
print.lp.bootgpd <- print.lp.gpd

################################################################################
## test.predict.gpd()

test.predict.gpd <- function(){
# no covariates
  u <- 14
  r.fit <- gpd(rain,th=u)
  co <- coef(r.fit)
  qgpd(0.9,exp(co[1]),co[2],u=u)

  checkEqualsNumeric(target=u,current = predict(r.fit,M=1/r.fit$rate)[[1]],msg="predict.gpd: retrieve threshold")

  checkEquals(target=predict(r.fit), current=rl(r.fit),msg="predict.gpd: predict with type=rl gives same as direct call to rl with default arguments")
  checkEquals(target=predict(r.fit,type="lp"), current=linearPredictors(r.fit),msg="predict.gpd: predict with type=rl gives same as direct call to rl with default arguments")

  t.fit <- r.fit
  t.fit$rate <- 1
  prob <- c(0.5,0.9,0.95,0.99,0.999)
  for(p in prob){
    checkEqualsNumeric(target = qgpd(p,exp(co[1]),co[2],u=u),
                       current = unlist(predict(t.fit,M=1/(1-p))),msg="predict.gpd: est ret levels no covariates")
  }

# with covariates

  n <- 1000
  M <- 1000

# boundary cases with single covariates in only one parameter

  xi <- 0.05
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,exp(X[,1]),xi)
  X$Y <- Y
  fit <- gpd(Y,data=X,phi=~a,th=0)
  co <- coef(fit)
  xi <- co[3]
  sig <- exp(cbind(rep(1,n),X[,1]) %*% co[1:2])
  pred <- predict(fit,M=M)

  checkEqualsNumeric(target=X$a,current=pred[[1]][,-1],msg="predict.gpd: ret level correct reporting of covariates for single cov in phi only")
  checkEqualsNumeric(target=qgpd(1-1/M,sig,xi,u=0),
                     current = pred[[1]][,1],msg="predict.gpd: ret level estimation with covariates in phi only")

  sig <- 2
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,sig,X[,2])
  X$Y <- Y
  fit <- gpd(Y,data=X,xi=~b,th=0)
  co <- coef(fit)
  sig <- exp(co[1])
  xi <- cbind(rep(1,n),X[,2]) %*% co[2:3]
  pred <- predict(fit,M=M)

  checkEqualsNumeric(target=X$b,current=pred[[1]][,-1],msg="predict.gpd: ret level correct reporting of covariates for single cov in xi only")
  checkEqualsNumeric(target=qgpd(1-1/M,sig,xi,u=0),
                     current = pred[[1]][,1],msg="predict.gpd: ret level estimation with covariates in xi only")

# covariates in both parameters

  n <- 1000
  M <- 1000
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,exp(X[,1]),X[,2])
  X$Y <- Y
  fit <- gpd(Y,data=X,phi=~a,xi=~b,th=0)
  co <- coef(fit)
  sig <- exp(cbind(rep(1,n),X[,1]) %*% co[1:2])
  xi <- cbind(rep(1,n),X[,2]) %*% co[3:4]

  checkEqualsNumeric(target=qgpd(1-1/M,sig,xi,u=0),
                     current = predict(fit,M=M)[[1]][,1],msg="predict.gpd: ret level estimation with covariates in phi and xi")

# check multiple M
  M <- c(10,50,100,500,1000)

  target <- sapply(M,function(m,sig,xi,u) qgpd(1-1/m,sig,xi,u=u),sig=sig,xi=xi,u=0)
  current <- predict(fit,M=M)

  for(i in 1:length(M)){
    checkEqualsNumeric(target[,i],current[[i]][,1],msg="predict.gpd: ret level estimation multiple M")
  }

# new data
  nx <- 20
  M <- 1000
  newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))

  sig <- exp(co[1] + newX[[1]] * co[2])
  xi <- co[3] + newX[[2]] * co[4]

  checkEqualsNumeric(target=qgpd(1-1/M,sigma=sig,xi=xi,u=0),current=predict(fit,M=M,newdata=newX)[[1]][,1],msg="predict.gpd: ret level ests with new data")

  checkEqualsNumeric(c(nx,5),dim(predict(fit,ci=TRUE,newdata=newX)[[1]]), msg="predict.gpd: dimension or return object for ci calc")
  checkEqualsNumeric(c(nx,4),dim(predict(fit,se=TRUE,newdata=newX)[[1]]), msg="predict.gpd: dimension or return object for se calc")
  checkEqualsNumeric(c(nx,6),dim(predict(fit,se=TRUE,ci=TRUE,newdata=newX)[[1]]), msg="predict.gpd: dimension or return object for se and ci calc")

  checkEquals(c("RL","2.5%","97.5%","se.fit","a","b"), colnames(predict(fit,se=TRUE,ci=TRUE)[[1]]), msg="predict.gpd: colnames of return obejct for se and ci calc, default alpha")
  checkEquals(c("RL","5%","95%","se.fit","a","b"), colnames(predict(fit,se=TRUE,ci=TRUE,alpha=0.1)[[1]]), msg="predict.gpd: colnames of return obejct for se and ci calc, alpha=0.1")

# alpha
  alpha <- c(0.01,0.05,0.1,0.2,0.5,0.9,0.99)
  z <- matrix(qnorm(c(alpha/2,1-alpha/2)),ncol=2)

  for(a in 1:length(alpha)){
    p <- predict(fit,alpha=alpha[a],ci=TRUE,se=TRUE)[[1]]
    checkEquals(current = colnames(p)[2:3],target = c(paste(100*alpha[a]/2,"%",sep=""),paste(100*(1-alpha[a]/2),"%",sep="")),msg="predict.gpd: labelling of confidence intervals")
    checkEqualsNumeric(target = p[,1] + z[a,1]*p[,4],current = p[,2], msg="predict.gpd: ret level Conf Interval calc for different alpha")
    checkEqualsNumeric(target = p[,1] + z[a,2]*p[,4],current = p[,3], msg="predict.gpd: ret level Conf Interval calc for different alpha")
  }

# linear predictors

  checkEqualsNumeric(target = cbind(co[1] + newX[[1]] * co[2],
                                    co[3] + newX[[2]] * co[4]),
                     current = predict(fit,newdata=newX,type="lp")[,1:2], msg="predict.gpd: linear predictors")

  checkEqualsNumeric(target = c(nx,6),dim(predict(fit,newdata=newX,se=TRUE,type="lp")), msg="predict.gpd: dimension of return object, linear predictor, se calc")
  checkEqualsNumeric(target = c(nx,8),dim(predict(fit,newdata=newX,ci=TRUE,type="lp")), msg="predict.gpd: dimension of return object, linear predictor, ci calc")

  name <- c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi")
  checkEquals(target = name, current = colnames(predict(fit,newdata=newX,ci=TRUE,type="lp"))[1:6],msg="predict.gpd: colnames for linear predictor return object")
  checkEquals(target = name, current = colnames(predict(fit,newdata=newX,ci=TRUE,se=TRUE,type="lp"))[1:6],msg="predict.gpd: colnames for linear predictor return object")

# unique

  newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1)) # checks for duplicates in one and both covariates.
  U <- !duplicated(newX)
  checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp"),
                     target = predict(fit,newdata=newX,unique.=FALSE,type="lp")[U,], msg="predict.gpd: functioning of argument unique, for linear predictor")
  checkEqualsNumeric(current = predict(fit,newdata=newX)[[1]],
                     target =  predict(fit,newdata=newX,unique.=FALSE)[[1]][U,], msg="predict.gpd: functioning of argument unique, for return levels")

# check standard errors - this takes a while since using bootstrap

  M <- c(10,100,500,1000,2000)
  newX <- data.frame("a"=c(1,2),"b"=c(-0.1,0.1))
  fit.p <- predict(fit, newdata=newX,se=TRUE,M=M)
  fit.seest <- unlist(lapply(fit.p,function(x) x[,2]))

  o <- options(warn=-1)
  fit.b <- bootgpd(fit,R=1000, trace=1100)
  options(o)
  fit.bp <- predict(fit.b,newdata=newX,all=TRUE,M=M)
  fit.seb <- lapply(fit.bp,function(X) apply(X,2,sd))
  fit.seboot <- unlist(fit.seb)

  checkEqualsNumeric(rep(0,length(fit.seboot)), (fit.seboot -  fit.seest) / fit.seest, tolerance=0.1,msg="predict.gpd: standard error estimate compared with bootstrap standard errors")
}

################################################################################
## test.predict.bgpd()

test.predict.bgpd <- function(){
# no covariates
  u <- 14
  r.fit <- gpd(rain,th=u,method="sim",trace=20000)

  checkEqualsNumeric(target=u,current=predict(r.fit,M=1/r.fit$map$rate)[[1]], msg="predict.bgpd: retrieve threshold")

  t.fit <- r.fit
  t.fit$map$rate <- 1
  p <- c(0.5,0.9,0.95,0.99,0.999)
  checkEqualsNumeric(target = sapply(p,function(p)mean(qgpd(p,exp(t.fit$param[,1]),t.fit$param[,2],u))),
                     current = unlist(predict(t.fit,M=1/(1-p))),msg="predict.bgpd: ret level estimation")

  checkEquals(target=predict(t.fit,M=1/(1-p)), current=rl(t.fit,M=1/(1-p)),msg="predict.bgpd: predict with type=rl gives same as direct call to rl with default arguments")
  checkEquals(target=predict(t.fit,type="lp"), current=linearPredictors(t.fit),msg="predict.bgpd: predict with type=rl gives same as direct call to rl with default arguments")
# with covariates

  n <- 100
  M <- 1000
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,exp(X[,1]),X[,2])
  X$Y <- Y
  fit <- gpd(Y,data=X,phi=~a,xi=~b,th=0,method="sim",trace=20000)

  sig <- apply(fit$param,1,function(v)exp(cbind(rep(1,n),X[,1]) %*% v[1:2]))
  xi <-  apply(fit$param,1,function(v)    cbind(rep(1,n),X[,2]) %*% v[3:4])

  qbgpd <- function(p,sig,xi,u)sapply(1:dim(sig)[1],function(i) mean(qgpd(p,sig[i,],xi[i,],u)))

  checkEqualsNumeric(target = qbgpd(1-1/M,sig,xi,0),
                     current = predict(fit,M=M)[[1]][,1],msg="predict.bgpd: ret level estimation with covariates")

# check multiple M
  M <- c(10,50,100,500,1000)

  temp1 <- sapply(M,function(m)qbgpd(1-1/m,sig,xi,u=0))
  temp2 <- predict(fit,M=M)

  for(i in 1:length(M)){
    checkEqualsNumeric(target = temp1[,i],current = temp2[[i]][,1],msg="predict.bgpd multiple M")
  }

# new data
  nx <- 20
  M <- 1000
  newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))

  sig <- exp(cbind(rep(1,nx),newX[,1]) %*% t(fit$param[,1:2]))
  xi <-      cbind(rep(1,nx),newX[,2]) %*% t(fit$param[,3:4])

  checkEqualsNumeric(target = as.matrix(newX), current = predict(fit,M=M,newdata=newX)[[1]][,2:3],msg="predict.bgpd : ret level estimation with new data, covariates aded correctly to output")
  checkEqualsNumeric(target = qbgpd(1-1/M,sig=sig,xi=xi,u=0),
                     current = predict(fit,M=M,newdata=newX)[[1]][,1],msg="predict.bgpd : ret level estimation with new data")

  p <- predict(fit,all=TRUE,newdata=newX)
  l <- lapply(p,function(l) apply(l,2,mean))
  m <- unlist(l)
  r <- predict(fit,newdata=newX)
  checkEqualsNumeric(target = m,current = r[[1]][,1],msg="predict.bgpd : ret level ests with new data")

  alpha <- c(0.05,0.1)
  for(a in alpha){
    l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
    l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
    m.L <- unlist(l.L)
    m.U <- unlist(l.U)
    r <- predict(fit,newdata=newX,ci=TRUE,alpha=a)
    checkEqualsNumeric(target = m.L,current = r[[1]][,3],msg="predict.bgpd : lower conf ints for ret levels with new data")
    checkEqualsNumeric(target = m.U,current = r[[1]][,4],msg="predict.bgpd : upper conf ints for ret levels with new data")
  }

# check linear predictors

  p <- predict(fit,type="lp",all=TRUE,newdata=newX)
  l <- lapply(p,function(l) apply(l,2,mean))
  m <- matrix(unlist(l),ncol=4,byrow=TRUE)
  r <- predict(fit,type="lp",newdata=newX)
  checkEqualsNumeric(target = m[,1:2],current = r[,1:2],msg="predict.bgpd : linear predictors of parameters with new data")

  alpha <- c(0.05,0.1)
  for(a in alpha){
    l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
    l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
    m.L <- matrix(unlist(l.L),ncol=4,byrow=TRUE)
    m.U <- matrix(unlist(l.U),ncol=4,byrow=TRUE)
    r <- predict(fit,type="lp",newdata=newX,ci=TRUE,alpha=a)
    checkEqualsNumeric(target = m.L[,1:2],current = r[,c(3,7)],msg="predict.bgpd : lower conf ints for linear predictors of parameters with new data")
    checkEqualsNumeric(target = m.U[,1:2],current = r[,c(4,8)],msg="predict.bgpd : upper conf ints for linear predictors of parameters with new data")
  }

  checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")[,1:2],
                     target = cbind(apply(cbind(rep(1,nx),newX[,1]) %*% t(fit$param[,1:2]),1,mean),
                                    apply(cbind(rep(1,nx),newX[,2]) %*% t(fit$param[,3:4]),1,mean)),
                     msg = "predict.bgpd: linear predictor estimates")

# structure of output

  checkEqualsNumeric(target = c(n,6), current = dim(predict(fit,ci=TRUE)[[1]]), msg="predict.bgpd: dimension of output with ci calculation")
  checkEqualsNumeric(target = c(n,6), current = dim(predict(fit,se=TRUE,ci=TRUE)[[1]]), msg="predict.bgpd: dimension of output with ci calculation")

  checkEquals(target = c("Mean", "50%","2.5%","97.5%","a","b"), colnames(predict(fit,ci=TRUE)[[1]]), msg="predict.bgpd: colnames of ret level ests with CI estimation")
  checkEquals(target = c("Mean", "50%","5%","95%","a","b"), colnames(predict(fit,alpha=0.1,ci=TRUE)[[1]]), msg="predict.bgpd: colnames of ret level ests with CI estimation, alpha=0.1")

  checkEqualsNumeric(target = c(nx,10), dim(predict(fit,newdata=newX,ci=TRUE,type="lp")), msg="predict.bgpd: dimension of linear predictor return object")

  cnames <- c("phi: Mean", "phi: 50%", "phi: 2.5%", "phi: 97.5%", "xi: Mean", "xi: 50%", "xi: 2.5%", "xi: 97.5%")# this specific format assumed by plot.rl.bgpd and plot.lp.bgpd

  checkEquals(current = cnames, target = colnames(predict(fit,newdata=newX,ci=TRUE,type="lp"))[1:8], msg="predict.bgpd: col names of lin predictors with CI calcs")
  checkEquals(current = cnames, target = colnames(predict(fit,newdata=newX,ci=TRUE,se=TRUE,type="lp"))[1:8], msg="predict.bgpd: col names of lin predictors with CI+SE calcs")

# unique
  newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1))

  checkEqualsNumeric(current = predict(fit,newdata=newX)[[1]], target = unique(predict(fit,newdata=newX,unique.=FALSE)[[1]]),msg="predict.bgpd: unique functioning for ret level ests")
  checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")[,], target = unique(predict(fit,newdata=newX,unique.=FALSE,type="lp")[,]),msg="predict.bgpd: unique functioning for lin pred ests")

}

################################################################################
## test.predict.bootgpd()

test.predict.bootgpd <- function(){
# functionality all tested already in test.predict.bgpd, so just check output of correct format.

  n <- 1000
  nx <- 9
  R <- 10
  nm <- 20

  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,exp(X[,1]),X[,2])
  X$Y <- Y
  fit <- gpd(Y,data=X,phi=~a,xi=~b,th=0)
  o <- options(warn=-1)
  boot <- bootgpd(fit,R=R,trace=100)
  options(o)

  newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
  from <- 10; to <- 500
  M <- seq(from,to,length=nm)
  pred <- predict(boot,newdata=newX,M=M,ci=TRUE)

  checkEquals(target=predict(boot), current=rl(boot),msg="predict.bootgpd: predict with type=rl gives same as direct call to rl with default arguments")
  checkEquals(target=predict(boot,type="lp"), current=linearPredictors(boot),msg="predict.bootgpd: predict with type=lp gives same as direct call to linearPredictors with default arguments")

  checkEqualsNumeric(target=nm,current=length(pred),msg="predict.bootgpd: output length")
  checkEquals(target=paste("M.",from,sep=""),current=names(pred)[1],msg="predict.bootgpd: names of output")
  checkEquals(target=paste("M.",to,sep=""),current=names(pred)[nm],msg="predict.bootgpd: names of output")

  cnames <- c( "Mean","50%","2.5%","97.5%",names(X)[1:2])
  checkEquals(target=cnames,current=colnames(pred[[1]]),msg="predict.bootgpd: colnames")

  checkEqualsNumeric(target=c(nx,6),current=dim(pred[[1]]),msg="predict.bootgpd: dimension")
  for(i in 1:nm){
    checkEqualsNumeric(target=newX[,1],current=pred[[i]][,5],msg="predict.bootgpd: covariates inoutput")
    checkEqualsNumeric(target=newX[,2],current=pred[[i]][,6],msg="predict.bootgpd: covariates inoutput")
  }

  par(mfrow=n2mfrow(nx))
  plot(pred,sameAxes=FALSE,type="median",main="Bootstrap median rl")
  plot(pred,sameAxes=FALSE,type="mean",main="Bootstrap mean rl")

}
