# Author: Harry Southworth
# Date: 2011-11-25
## Purpose: Create a predict method for objects of class evmOpt, evmSim
##          and evmBoot that
##          returns parameters, return levels or (maybe) return periods,
##          depending on arguments given.
#
# predict.evm
# predict.evmSim
# predict.evmBoot
# rl
# rl.evm
# rl.evmSim
# rl.evmBoot
# linearPredictors
# linearPredictors.evm
# linearPredictors.evmSim
# linearPredictors.evmBoot

################################################################################
## evm

predict.evmOpt <-
    # Get predictions for an evm object. These can either be the linear predictors
    # or return levels.
function(object, M=1000, newdata=NULL, type="return level", se.fit=FALSE,
         ci.fit=FALSE, alpha=.050, unique.=TRUE, ...){
    theCall <- match.call()

    res <- switch(type,
                  "rl"=, "return level" = rl.evmOpt(object, M, newdata,
                                                 se.fit=se.fit, ci.fit=ci.fit,
                                                 alpha=alpha, unique.=unique.),
                  "lp" =,"link" = linearPredictors.evmOpt(object, newdata, se.fit,
                                                   ci.fit, alpha, unique.=unique.)
                  )
    res
}

## Linear predictor functions for GPD

linearPredictors.evmOpt <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                             alpha=.050, unique.=TRUE, full.cov=FALSE, ...){

    D <- texmexMakeNewdataD(object, newdata)

    if (unique.){
        z <- do.call('cbind', D)
        u <- !duplicated(z)
        D <- lapply(D, function(x, u) {
                           if(is.matrix(x[u,]))  x[u, ]
                           else if(ncol(x) == 1)  cbind(x[u,])
                           else t(cbind(x[u,]))
                       }, u=u )
    }

    res <- texmexMakeParams(coef(object), D)
    colnames(res) <- names(D)
browser()
    # Get the covariance matrices - one for every unique observation
    if(ci.fit | se.fit | full.cov){
      cov.se <- texmexMakeCovariance(object$cov, D)
      # Get standard errors
      ses <- t(sapply(cov.se, function(x){ sqrt(diag(x)) }))
      colnames(ses) <- paste(colnames(res), '.se', sep = '')
    }

    if (ci.fit){
        ci <- texmexMakeCI(res, ses, alpha)
        res <- cbind(res, ci)
    } # Close if(ci.fit

    if (se.fit){
        res <- cbind(res, ses)
    } # Close if(se.fit

    for (i in 1:length(D)){
      res <- addCov(res, D[[i]])
    }

    res <- list(link=res,family=object$family)
    
    if (full.cov){
        res$cov <- cov.se
    }

    oldClass(res) <- "lp.evmOpt"
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


rl.evmOpt <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                       alpha=.050, unique.=TRUE, ...){
    co <- linearPredictors.evmOpt(object, newdata=newdata, unique.=unique., full.cov=TRUE)
    covs <- co$cov # list of covariance matrices, one for each (unique) observation
    co <- co$link
    X <- co[,-(1:length(object$data$D))]
    if(is.null(dim(X))){
      X <- matrix(X)
      dimnames(X) <- list(dimnames(co)[[1]],dimnames(co)[[2]][-(1:length(object$data$D))])
    }

    delta <- object$family$delta
    rl <- object$family$rl

    res <- lapply(M, rl, param=co, model=object)

    getse <- function(o, co, M, delta, covs){
        dxm <- lapply(split(co, 1:nrow(co)), delta, m=M, model=o)

        # Get (4.15) of Coles, page 82, adjusted for phi = log(sigma)
        se <- sapply(1:length(covs),
                     function(i, dxm, covs){
                        covs <- covs[[i]]; dxm <- c(dxm[[i]])
                        sqrt(mahalanobis(dxm, center=rep(0, ncol(covs)), cov=covs, inverted=TRUE))
                     }, dxm=dxm, covs=covs)
        se
    }

    if (ci.fit){
        ci.fun <- function(i, object, co, M, res, alpha, delta, covs){
            wh <- res[[i]];
            se <- getse(object, co, M[i], delta=delta, covs=covs)
            lo <- wh - qnorm(1 - alpha/2)*se
            hi <- wh + qnorm(1 - alpha/2)*se
            wh <- cbind(wh, lo=lo, hi=hi)

            colnames(wh) <- c("RL", paste(100*alpha/2, "%", sep = ""),
                              paste(100*(1 - alpha/2), "%", sep = ""))
            wh
        } # ci.fun
        res <- lapply(1:length(M), ci.fun, object=object, co=co,
                                           M=M, res=res, alpha=alpha,
                                           delta=delta, covs=covs)
    } # Close if (ci.fit

    if (se.fit){
        se.fun <- function(i, object, co, M, res, alpha, delta, covs){
            wh <- res[[i]]
            se <- getse(object, co, M[i], delta=delta, covs=covs)
            wh <- cbind(RL=wh, se.fit=se)
            wh
        } # ci.fun
        res <- lapply(1:length(M), se.fun, object=object, co=co,
                                           M=M, res=res, alpha=alpha,
                                           delta=delta, covs=covs)
    }

    cov.fun <- function(i,res){
      wh <- res[[i]]
      wh <- addCov(wh,X)
      wh
    }
    res <- lapply(1:length(M), cov.fun,res=res)

    names(res) <- paste("M.", M, sep = "")
    oldClass(res) <- "rl.evmOpt"
    res
}

################################################################################
## evmSim

predict.evmSim <- function(object, M=1000, newdata=NULL, type="return level",
                         se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE,
                         all=FALSE, sumfun=NULL, ...){
    theCall <- match.call()

    res <- switch(type,
                  "rl" = , "return level" = rl.evmSim(object, M=M, newdata=newdata,
                                                    se.fit=se.fit, ci.fit=ci.fit,
                                                    alpha=alpha, unique.=unique., all=all,
                                                    sumfun=sumfun,...),
                  "lp" = , "link" = linearPredictors.evmSim(object, newdata=newdata,
                                                      se.fit=se.fit, ci.fit=ci.fit,
                                                      alpha=alpha, unique.=unique., all=all,
                                                      sumfun=sumfun,...)
                  )
    res
}

linearPredictors.evmSim <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                                     alpha=.050, unique.=TRUE, all=FALSE, sumfun=NULL, ...){
    if (se.fit){ warning("se.fit not implemented - ignoring") }

    D <- texmexMakeNewdataD(object$map, newdata)

    X.all <- do.call("cbind", D)
    ModelHasCovs <- ncol(X.all) > length(D)

    if(ModelHasCovs){
      covCols <- apply(X.all, 2, function(x) !all(x==1))
      Xnames <- colnames(X.all)
      if(sum(covCols) == 1){
        X.all <- X.all[, covCols, drop=FALSE]
      }
      else {
        X.all <- X.all[, covCols]
      }
      colnames(X.all) <- Xnames[covCols]
    }

    if (unique.){
        X.all <- unique(X.all)
        u <- !duplicated(do.call("cbind", D))
        D <- lapply(D, function(x, u){ x[u,, drop=FALSE] }, u=u)
    }

    # Get matrices of parameters (i.e. split full parameter matrix into phi, xi whatever)
    param <- texmexGetParam(D, object$param)

    # Get linear predictors
    res <- lapply(1:nrow(D[[1]]), # For each observation get matrix of parameters
              function(i, x, p){
                  wh <- lapply(1:length(D),
                               function(j, x, p, i){
                                   rowSums(t(t(p[[j]]) * c(x[[j]][i, ])))
                               }, x=x, p=p, i=i)
                  wh <- do.call("cbind", wh)
                  colnames(wh) <- names(x)
                  wh
                }, x=D, p=param)
    # res should be a list containing a matrix for each observation.
    # The matrix represents the simulated posterior, one column for each
    # major parameter (i.e. linear predictors)

    ############################################################################
    ## Hard part should be done now. Just need to summarize

    if (ci.fit){
        # Need to get names by pasting together CI names and parameter names
        wh <- texmexMakeCISim(res[[1]], alpha=alpha, object=object, sumfun=sumfun)
        wh <- colnames(wh)
        wh <- paste(rep(names(D), ea=length(wh)), wh, sep = ": ")

        # Need to get order of output correct, so need to faff about with transposing
        res <- sapply(res, function(x){
                               t(texmexMakeCISim(x, alpha=alpha, object=object, sumfun=sumfun))
                           })
        res <- t(res)
        colnames(res) <- wh
    }

    else if (all){ res <- res }
    else { # Just point estimates
        res <- t(sapply(res, function(x){ apply(x, 2, mean) }))
    }
    if(!all){
      if(ModelHasCovs){
        for (i in 1:length(D)){
            res <- addCov(res,D[[i]])
        }
      }
    }
    else {
        if (ModelHasCovs & nrow(X.all) != length(res)){
            stop("Number of unique combinations of covariates doesn't match the number of parameters")
        }
        for (i in 1:length(res)){
          if(ModelHasCovs){
            res[[i]] <- cbind(res[[i]], matrix(rep(X.all[i,], nrow(res[[i]])),
                                               nrow=nrow(res[[i]]), byrow=TRUE))
            colnames(res[[i]]) <- c(names(D), colnames(X.all))
          } else {
            colnames(res[[i]]) <- names(D)
          }
        }
    }

    res <- list(link=res,family=object$map$family)
    oldClass(res) <- "lp.evmSim"
    res
}

rl.evmSim <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE, all=FALSE, sumfun=NULL,...){
    if (se.fit){ warning("se.fit not implemented") }

    co <- linearPredictors.evmSim(object, newdata=newdata, unique.=unique., all=TRUE, sumfun=NULL)$link
    # XXX Next line seems silly! Why not compute it from the line above?
    Covs <- linearPredictors.evmSim(object, newdata=newdata, unique.=unique., sumfun=NULL)$link
    X <- Covs[,-(1:length(object$map$data$D))]
    if(is.null(dim(X))){
      X <- matrix(X)
      dimnames(X) <- list(dimnames(Covs)[[1]],dimnames(Covs)[[2]][-(1:length(object$map$data$D))])
    }

    sim.rl <- function(m, param, model){
        rl <- model$family$rl
        rl(m=m, param, model)
    }

    # co is a list with one element for each unique item in
    # new data. Need to loop over vector M and the elements of co

    getrl <- function(m, co, ci.fit, alpha, all, object){
        res <- sapply(co, sim.rl, m=m, model=object$map)
        if (ci.fit){
            res <- texmexMakeCISim(res, alpha, object$map, sumfun, M=m)
        } # Close if (ci.fit
        else if (!all){
            res <- apply(res, 2, mean)
        }
        res
    }

    res <- lapply(M, getrl, co=co, ci.fit=ci.fit, alpha=alpha, all=all, object=object)

    if(!all){
      cov.fun <- function(i,res){
        wh <- res[[i]]
        wh <- addCov(wh,X)
        wh
      }
      res <- lapply(1:length(M), cov.fun,res=res)
    }

    names(res) <- paste("M.", M, sep = "")
    oldClass(res) <- "rl.evmSim"
    res
}

################################################################################
## evmBoot

predict.evmBoot <- function(object, M=1000, newdata=NULL, type="return level",
                            se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE,
                            all=FALSE, sumfun=NULL, ...){
    theCall <- match.call()

    res <- switch(type,
                  "rl" = , "return level" = rl.evmBoot(object, newdata=newdata, M=M,
                                                       se.fit=se.fit, ci.fit=ci.fit,
                                                       alpha=alpha, unique.=TRUE,
                                                       all=all, sumfun=sumfun,...),
                  "lp" = , "link" = linearPredictors.evmBoot(object, newdata=newdata,
                                                         se.fit=se.fit, ci.fit=ci.fit,
                                                         alpha=alpha, unique.=TRUE,
                                                         all=all, sumfun=sumfun,...)
                  )
    res
}

namesBoot2sim <- function(bootobject){
    names(bootobject) <- c("call", "param", "map")
    bootobject
}

linearPredictors.evmBoot <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=.050,
                                 unique.=TRUE, all=FALSE, sumfun=NULL,...){
    # This should just be the same as for an evmSim object, but some
    # names and stuff are different.
  object <- namesBoot2sim(object)
  res <- linearPredictors.evmSim(object, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit, all=all, unique.=unique., alpha=alpha, sumfun=sumfun,...)
  oldClass(res) <- "lp.evmBoot"
  res
}

rl.evmBoot <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, alpha=0.050, unique.=TRUE, all=FALSE, sumfun=NULL,...){
    # This should just be the same as for an evmSim object, but some
    # names are different.
  object <- namesBoot2sim(object)
  res <- rl.evmSim(object, M=M, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit,alpha=alpha, unique.=unique., all=all, sumfun=sumfun,...)
  oldClass(res) <- "rl.evmBoot"
  res
}

################################################################################
## Method functions

print.rl.evmOpt <- function(x, digits=3, ...){
    nms <- names(x)
    newnms <- paste("M =", substring(nms, 3), "predicted return level:\n")
    lapply(1:length(x), function(i, o, title){
                                 cat(title[i])
                                 print(o[[i]], digits=digits,...)
                                 cat("\n")
                                 NULL}, o=x, title=newnms)
    invisible(x)
}

summary.rl.evmOpt <- function(object, digits=3, ...){
    print.rl.evmOpt(object, digits=digits, ...)
}

print.rl.evmSim    <- print.rl.evmOpt
print.rl.evmBoot <- print.rl.evmOpt

summary.rl.evmSim    <- summary.rl.evmOpt
summary.rl.evmBoot <- summary.rl.evmOpt


print.lp.evmOpt <- function(x, digits=3, ...){
    cat("Linear predictors:\n")
    print(unclass(x$link), digits=3,...)
    invisible(x)
}

summary.lp.evmOpt <- function(object, digits=3, ...){
    print.lp.evmOpt(object, digits=3, ...)
}

summary.lp.evmSim    <- summary.lp.evmOpt
summary.lp.evmBoot <- summary.lp.evmOpt

print.lp.evmSim    <- print.lp.evmOpt
print.lp.evmBoot <- print.lp.evmOpt

################################################################################
## test.predict.evm()

test.predict.evmOpt <- function(){

  for(Family in list(gpd,gev)){
    set.seed(20130513)
    pst <- function(msg) texmexPst(msg,Family=Family)
    
# no covariates

    u    <- switch(Family$name,GPD=14,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    r.fit <- evm(data,th=u,family=Family)
    co <- coef(r.fit)

    if(Family$name == "GPD")checkEqualsNumeric(target=u,current = predict(r.fit,M=1/r.fit$rate)[[1]],msg=pst("predict.evmOpt: GPD retrieve threshold"))

    checkEquals(target=predict(r.fit), current=rl(r.fit),msg=pst("predict.evmOpt: predict with type=rl gives same as direct call to rl with default arguments"))
    checkEquals(target=predict(r.fit,type="lp"), current=linearPredictors(r.fit),msg=pst("predict.evmOpt: predict with type=lp gives same as direct call to linearpredictors with default arguments"))

    r.fit$rate <- 1
    prob <- c(0.5,0.9,0.95,0.99,0.999)
    for(p in prob){
      checkEqualsNumeric(target = Family$quant(p,t(co),r.fit),
                         current = unlist(predict(r.fit,M=1/(1-p))),msg=pst("predict.evmOpt: est ret levels no covariates"))
    }
  
# with covariates

    n <- 1000
    M <- 1000

# boundary cases with single covariates in only one parameter
    mu <- 1
    xi <- 0.05
    
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=cbind(X[,1],xi),GEV=cbind(mu,X[,1],xi))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,phi=~a,th=th,family=Family)
    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]
    pred <- predict(fit,M=M)

    checkEqualsNumeric(target=X$a,current=pred[[1]][,-1],msg=pst("predict.evmOpt: ret level correct reporting of covariates for single cov in phi only"))
    checkEqualsNumeric(target=Family$quant(1-1/M,AllCo,fit),
                       current = pred[[1]][,1],msg=pst("predict.evmOpt: ret level estimation with covariates in phi only"))

    mu <- 1
    sig <- 2

    param <- switch(Family$name,GPD=cbind(log(sig),X[,2]),GEV=cbind(mu,log(sig),X[,2]))
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,xi=~b,th=th,family=Family)
    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]
    pred <- predict(fit,M=M)

    checkEqualsNumeric(target=X$b,current=pred[[1]][,-1],msg=pst("predict.evmOpt: ret level correct reporting of covariates for single cov in xi only"))
    checkEqualsNumeric(target=Family$quant(1-1/M,AllCo,fit),
                       current = pred[[1]][,1],msg=pst("predict.evmOpt: ret level estimation with covariates in xi only"))

# covariates in all parameters

    param <- switch(Family$name,GPD=cbind(X[,1],X[,2]),GEV=cbind(X[,1],X[,1],X[,2]))
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- switch(Family$name,
                  GPD = evm(Y, data=X,        phi=~a, xi=~b, th=th,family=Family),
                  GEV = evm(Y, data=X, mu=~a, phi=~a, xi=~b, th=th,family=Family))  

    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]

    checkEqualsNumeric(target=Family$quant(1-1/M,AllCo,fit),
                      current = predict(fit,M=M)[[1]][,1],msg=pst("predict.evmOpt: ret level estimation with covariates in all parameters"))

# check multiple M
    M <- c(10,50,100,500,1000)

    target <- sapply(M,function(m,AllCo,fit) Family$quant(1-1/m,AllCo,fit),AllCo=AllCo,fit=fit)
    current <- predict(fit,M=M)

    for(i in 1:length(M)){
      checkEqualsNumeric(target[,i],current[[i]][,1],msg=pst("predict.evmOpt: ret level estimation multiple M"))
    }

# new data
    nx <- 20
    M <- 1000
    newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
    AllCoNew <- predict(fit,type="lp",newdata=newX)$link[,1:length(Family$param)]

    checkEqualsNumeric(target=Family$quant(1-1/M,AllCoNew,fit),current=predict(fit,M=M,newdata=newX)[[1]][,1],msg=pst("predict.evmOpt: ret level ests with new data"))

    checkEqualsNumeric(dim(AllCoNew) + c(0,3),dim(predict(fit,ci=TRUE,newdata=newX)[[1]]), msg=pst("predict.evmOpt: dimension of return object for ci calc"))
    checkEqualsNumeric(dim(AllCoNew) + c(0,2),dim(predict(fit,se=TRUE,newdata=newX)[[1]]), msg=pst("predict.evmOpt: dimension of return object for se calc"))
    checkEqualsNumeric(dim(AllCoNew) + c(0,4),dim(predict(fit,se=TRUE,ci=TRUE,newdata=newX)[[1]]), msg=pst("predict.evmOpt: dimension of return object for se and ci calc"))

    Labels <- switch(Family$name,GPD=c("a","b"),GEV=c("a","a","b"))
    
    checkEquals(c("RL","2.5%","97.5%","se.fit",Labels), colnames(predict(fit,se=TRUE,ci=TRUE)[[1]]), msg=pst("predict.evmOpt: colnames of return obejct for se and ci calc, default alpha"))
    checkEquals(c("RL","5%","95%","se.fit",Labels), colnames(predict(fit,se=TRUE,ci=TRUE,alpha=0.1)[[1]]), msg=pst("predict.evmOpt: colnames of return obejct for se and ci calc, alpha=0.1"))

# alpha

    alpha <- c(0.01,0.05,0.1,0.2,0.5,0.9,0.99)
    z <- matrix(qnorm(c(alpha/2,1-alpha/2)),ncol=2)

    for(a in 1:length(alpha)){
      p <- predict(fit,alpha=alpha[a],ci=TRUE,se=TRUE)[[1]]
      checkEquals(current = colnames(p)[2:3],target = c(paste(100*alpha[a]/2,"%",sep=""),paste(100*(1-alpha[a]/2),"%",sep="")),msg=pst("predict.evmOpt: labelling of confidence intervals"))
      checkEqualsNumeric(target = p[,1] + z[a,1]*p[,4],current = p[,2], msg=pst("predict.evmOpt: ret level Conf Interval calc for different alpha"))
      checkEqualsNumeric(target = p[,1] + z[a,2]*p[,4],current = p[,3], msg=pst("predict.evmOpt: ret level Conf Interval calc for different alpha"))
    }


# linear predictors

    checkEqualsNumeric(target = c(nx,switch(Family$name,GPD=6,GEV=9)), dim(predict(fit,newdata=newX,se=TRUE,type="lp")$link), msg=pst("predict.evmOpt: dimension of return object, linear predictor, se calc"))
    checkEqualsNumeric(target = c(nx,switch(Family$name,GPD=8,GEV=12)),dim(predict(fit,newdata=newX,ci=TRUE,type="lp")$link), msg=pst("predict.evmOpt: dimension of return object, linear predictor, ci calc"))

    nameCI <- switch(Family$name,GPD = c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi","a","b"),GEV = c("mu","phi", "xi", "mu.lo","mu.hi","phi.lo", "phi.hi", "xi.lo", "xi.hi","a","a","b"))
    nameSE <- switch(Family$name,GPD = c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi","phi.se", "xi.se","a","b"),GEV = c("mu","phi", "xi", "mu.lo","mu.hi","phi.lo", "phi.hi", "xi.lo", "xi.hi","mu.se","phi.se", "xi.se","a","a","b"))
    checkEquals(target = nameCI, current = colnames(predict(fit,newdata=newX,ci=TRUE,type="lp")$link),msg=pst("predict.evmOpt: colnames for linear predictor return object"))
    checkEquals(target = nameSE, current = colnames(predict(fit,newdata=newX,ci=TRUE,se=TRUE,type="lp")$link),msg=pst("predict.evmOpt: colnames for linear predictor return object"))

# unique

    newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1)) # checks for duplicates in one and both covariates.
    U <- !duplicated(newX)
    checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")$link,
                       target = predict(fit,newdata=newX,unique.=FALSE,type="lp")$link[U,], msg=pst("predict.evmOpt: functioning of argument unique, for linear predictor"))
    checkEqualsNumeric(current = predict(fit,newdata=newX)[[1]],
                       target =  predict(fit,newdata=newX,unique.=FALSE)[[1]][U,], msg=pst("predict.evmOpt: functioning of argument unique, for return levels"))
    
# check standard errors - this takes a while since using bootstrap

    M <- c(10,100,500,1000,2000)
    newX <- data.frame("a"=rep(c(1,-1,2,-2),2),"b"=c(rep(0.1,4),rep(-0.1,4)))
    fit.p <- predict(fit, newdata=newX,se=TRUE,M=M)
    fit.seest <- unlist(lapply(fit.p,function(x) x[,2]))
  
    o <- options(warn=-1)
    fit.b <- evmBoot(fit,R=1000, trace=1100)
    options(o)
    fit.bp <- predict(fit.b,newdata=newX,all=TRUE,M=M)
    fit.seb <- lapply(fit.bp,function(X) apply(X,2,sd))
    fit.seboot <- unlist(fit.seb)

    checkTrue(all(abs((fit.seboot -  fit.seest) / fit.seest) < 0.3),msg=pst("predict.evmOpt: return level standard error estimate compared with bootstrap standard errors"))
  }  
}

################################################################################
## test.predict.evmSim()

test.predict.evmSim <- function(){

  for(Family in list(gpd,gev)){
    
    pst <- function(msg) texmexPst(msg,Family=Family)
    set.seed(20130513)
    # no covariates
    
    u    <- switch(Family$name,GPD=14,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)

    r.fit <- evm(data,th=u,family=Family,method="sim",trace=50000)
    co <- coef(r.fit)
   
    if(Family$name == "GPD"){
      checkEqualsNumeric(target=u,current=predict(r.fit,M=1/r.fit$map$rate)[[1]], msg=pst("predict.evmSim: retrieve threshold"))
    }

    r.fit$map$rate <- 1
    p <- c(0.5,0.9,0.95,0.99,0.999)
    checkEqualsNumeric(target = Family$quant(p,t(co),r.fit$map), tolerance=0.01,
                       current = unlist(predict(r.fit,M=1/(1-p))),msg=pst("predict.evmSim: ret level estimation"))

    checkEquals(target=predict(r.fit,M=1/(1-p)), current=rl(r.fit,M=1/(1-p)),msg=pst("predict.evmSim: predict with type=rl gives same as direct call to rl with default arguments"))
    checkEquals(target=predict(r.fit,type="lp")$link, current=linearPredictors(r.fit)$link,msg=pst("predict.evmSim: predict with type=rl gives same as direct call to rl with default arguments"))

# with covariates

    n <- 1000
    M <- 1000
    mu <- 1
  
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=X,GEV=cbind(mu,X))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    start <- switch(Family$name,GPD=c(0,1,0,1),GEV=c(1,0,1,0,1))
    fit <- evm(Y,data=X,phi=~a,xi=~b,th=th,method="sim",trace=50000,family=Family,start=start)

    AllCo <- predict(fit,type="lp",all=TRUE)$link
    PostMeanRL <- function(AllCo,M){
      AllQuant <- lapply(AllCo,function(X)Family$quant(1-1/M,X[,1:length(Family$param)],fit$map))
      sapply(AllQuant,mean)
    }
 
    checkEqualsNumeric(target=PostMeanRL(AllCo,M),current = predict(fit,M=M)[[1]][,1],msg=pst("predict.evmSim: ret level estimation with covariates in all parameters"))
# multiple M
  
    M <- c(10,50,100,500,1000)
  
    target <- lapply(M,function(m) PostMeanRL(AllCo,m))
    current <- predict(fit,M=M)
  
    for(i in 1:length(M)){
      checkEqualsNumeric(target[[i]],current[[i]][,1],tolerance=0.02,msg=pst("predict.evmSim: ret level estimation multiple M"))
    }

# new data
    nx <- 20
    M <- 1000
    newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))

    AllCoNew <- predict(fit,type="lp",newdata=newX,all=TRUE)$link
    
    checkEqualsNumeric(target=PostMeanRL(AllCoNew,M),current=predict(fit,M=M,newdata=newX)[[1]][,1],msg=pst("predict.evmSim: ret level ests with new data"))
    checkEqualsNumeric(target = as.matrix(newX), current = predict(fit,M=M,newdata=newX)[[1]][,2:3],msg=pst("predict.evmSim: ret level estimation with new data, covariates added correctly to output"))


    p <- predict(fit,all=TRUE,newdata=newX)
    alpha <- c(0.05,0.1)
    for(a in alpha){
      l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
      l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
      m.L <- unlist(l.L)
      m.U <- unlist(l.U)
      r <- predict(fit,newdata=newX,ci=TRUE,alpha=a)
      checkEqualsNumeric(target = m.L,current = r[[1]][,3],msg=pst("predict.evmSim : lower conf ints for ret levels with new data"))
      checkEqualsNumeric(target = m.U,current = r[[1]][,4],msg=pst("predict.evmSim : upper conf ints for ret levels with new data"))
    }

# check linear predictors
    p <- predict(fit,type="lp",all=TRUE,newdata=newX)$link
    l <- lapply(p,function(l) apply(l,2,mean))
    m <- matrix(unlist(l),ncol=length(l[[1]]),byrow=TRUE)
    r <- predict(fit,type="lp",newdata=newX)$link
    checkEqualsNumeric(target = m,current = r,msg=pst("predict.evmSim : linear predictors of parameters with new data"))
    Offset <- switch(Family$name,GEV=1,GPD=0)
    checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")$link[,Offset +1:2],
                       target = cbind(apply(cbind(rep(1,nx),newX[,1]) %*% t(fit$param[,Offset + 1:2]),1,mean),
                                      apply(cbind(rep(1,nx),newX[,2]) %*% t(fit$param[,Offset + 3:4]),1,mean)),
                       msg = pst("predict.evmSim: linear predictor estimates"))
    
    alpha <- c(0.05,0.1)
    for(a in alpha){
      l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
      l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
      m.L <- matrix(unlist(l.L),ncol=length(l[[1]]),byrow=TRUE)
      m.U <- matrix(unlist(l.U),ncol=length(l[[1]]),byrow=TRUE)
      r <- predict(fit,type="lp",newdata=newX,ci=TRUE,alpha=a)$link
      npar <- length(Family$param)
      checkEqualsNumeric(target = m.L[,1:npar],current = r[,(1:(4*npar))[rep(c(F,F,T,F),npar)]],msg=pst("predict.evmSim : lower conf ints for linear predictors of parameters with new data"))
      checkEqualsNumeric(target = m.U[,1:npar],current = r[,(1:(4*npar))[rep(c(F,F,F,T),npar)]],msg=pst("predict.evmSim : upper conf ints for linear predictors of parameters with new data"))
    }

# structure of output

    checkEqualsNumeric(target = c(n,6), current = dim(predict(fit,ci=TRUE)[[1]]), msg=pst("predict.evmSim: dimension of output with ci calculation"))
    o <- options(warn=-1) # since se=TRUE gives a warning
    checkEqualsNumeric(target = c(n,6), current = dim(predict(fit,se=TRUE,ci=TRUE)[[1]]), msg=pst("predict.evmSim: dimension of output with ci and se calculation"))
    options(o)
    
    checkEquals(target = c("Mean", "50%","2.5%","97.5%","a","b"), colnames(predict(fit,ci=TRUE)[[1]]), msg=pst("predict.evmSim: colnames of ret level ests with CI estimation"))
    checkEquals(target = c("Mean", "50%","5%","95%","a","b"), colnames(predict(fit,alpha=0.1,ci=TRUE)[[1]]), msg=pst("predict.evmSim: colnames of ret level ests with CI estimation, alpha=0.1"))

    checkEqualsNumeric(target = c(nx,4*npar+2), dim(predict(fit,newdata=newX,ci=TRUE,type="lp")$link), msg=pst("predict.evmSim: dimension of linear predictor return object"))

    cnamesGPD <- c("phi: Mean", "phi: 50%", "phi: 2.5%", "phi: 97.5%", "xi: Mean", "xi: 50%", "xi: 2.5%", "xi: 97.5%")# this specific format assumed by plot.rl.evmSim and plot.lp.evmSim
    cnamesGEV <- c("mu: Mean",  "mu: 50%",  "mu: 2.5%",  "mu: 97.5%",  "phi: Mean", "phi: 50%", "phi: 2.5%", "phi: 97.5%", "xi: Mean", "xi: 50%", "xi: 2.5%", "xi: 97.5%")
    cnames<- switch(Family$name,GPD=cnamesGPD,GEV=cnamesGEV)
  
    checkEquals(current = cnames, target = colnames(predict(fit,newdata=newX,ci=TRUE,type="lp")$link)[1:(4*npar)], msg=pst("predict.evmSim: col names of lin predictors with CI calcs"))
    o <- options(warn=-1) # since se=TRUE gives a warning
    checkEquals(current = cnames, target = colnames(predict(fit,newdata=newX,ci=TRUE,se=TRUE,type="lp")$link)[1:(4*npar)], msg=pst("predict.evmSim: col names of lin predictors with CI+SE calcs"))
    options(o)
    
# unique
    newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1))

    checkEqualsNumeric(current = predict(fit,newdata=newX)[[1]], target = unique(predict(fit,newdata=newX,unique.=FALSE)[[1]]),pst("predict.evmSim: unique functioning for ret level ests"))
    checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")$link[,], target = unique(predict(fit,newdata=newX,unique.=FALSE,type="lp")$link[,]),msg=pst("predict.evmSim: unique functioning for lin pred ests"))
  }
}

################################################################################
## test.predict.evmBoot()

test.predict.evmBoot <- function(){
# functionality all tested already in test.predict.evm, so just check output of correct format.

  n <- 1000
  nx <- 9
  R <- 10
  nm <- 20

  for(Family in list(gpd,gev)){
    
    pst <- function(msg) texmexPst(msg,Family=Family)
    set.seed(20130513)
    
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=X,GEV=cbind(5,X))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,phi=~a,xi=~b,th=th,family=Family)

    o <- options(warn=-1)
    boot <- evmBoot(fit,R=R,trace=100)
    options(o)

    newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
    from <- 10; to <- 500
    M <- seq(from,to,length=nm)
    pred <- predict(boot,newdata=newX,M=M,ci=TRUE)

    checkEquals(target=predict(boot), current=rl(boot),msg=pst("predict.evmBoot: predict with type=rl gives same as direct call to rl with default arguments"))
    checkEquals(target=predict(boot,type="lp"), current=linearPredictors(boot),msg=pst("predict.evmBoot: predict with type=lp gives same as direct call to linearPredictors with default arguments"))

    checkEqualsNumeric(target=nm,current=length(pred),msg=pst("predict.evmBoot: output length"))
    checkEquals(target=paste("M.",from,sep=""),current=names(pred)[1],msg=pst("predict.evmBoot: names of output"))
    checkEquals(target=paste("M.",to,sep=""),current=names(pred)[nm],msg=pst("predict.evmBoot: names of output"))

    cnames <- c( "Mean","50%","2.5%","97.5%",names(X)[1:2])
    checkEquals(target=cnames,current=colnames(pred[[1]]),msg=pst("predict.evmBoot: colnames"))

    checkEqualsNumeric(target=c(nx,6),current=dim(pred[[1]]),msg=pst("predict.evmBoot: dimension"))
    for(i in 1:nm){
      checkEqualsNumeric(target=newX[,1],current=pred[[i]][,5],msg=pst("predict.evmBoot: covariates inoutput"))
      checkEqualsNumeric(target=newX[,2],current=pred[[i]][,6],msg=pst("predict.evmBoot: covariates inoutput"))
    }

    par(mfrow=n2mfrow(nx))
    plot(pred,sameAxes=FALSE,type="median",main="Bootstrap median rl")
    plot(pred,sameAxes=FALSE,type="mean",main="Bootstrap mean rl")
  }
}
