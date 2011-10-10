bootmex <- 
    # Bootstrap inference for a conditional multivaratiate extremes model.
function (x, R = 100, nPass = 3, trace = 10) {
    theCall <- match.call()
    if (class(x) != "mex"){
      stop("object must be of type 'mex'")
    }
    
# Construct the object to be returned.
    ans <- list()
    ans$call <- theCall
    
    getTran <- function(i, x, data, mod, th, qu, margins) {
        x <- c(x[, i])
        param <- mod[[i]]$coefficients
        th <- th[i]
        qu <- qu[i]
        data <- c(data[, i])
        if (margins == "gumbel"){
			     res <- revTransform(exp(-exp(-x)), data = data, th = th, 
            					qu = qu, sigma = exp(param[1]), xi = param[2])
		    } else {
			    y <- x
			    y[y < 0] <- exp(y[y < 0]) / 2
			    y[y >= 0] <- 1 - exp(-y[y >= 0]) / 2
			    res <- revTransform(y, data = data, th = th, 
					     			qu = qu, sigma = exp(param[1]), xi = param[2])
		    }
        res
    }

    mar <- x$margins
    dep <- x$dependence
    which <- dep$which
    constrain <- dep$constrain
    v <- dep$v
    dqu <- dep$dqu
    dth <- dep$dth
    margins <- dep$margins        
    penalty <- mar$penalty
    priorParameters <- mar$priorParameters

    n <- dim(mar$transformed)[[1]]
    d <- dim(mar$transformed)[[2]]
    dqu <- rep(dqu, d) 
    dependent <- (1:d)[-which]
    
    ans$simpleDep <- dep$coefficients
    ans$dqu <- dqu
    ans$which <- which
    ans$R <- R
    ans$simpleMar <- mar
    ans$margins <- margins
    ans$constrain <- constrain

    innerFun <- function(i, x, which, dth, dqu, margins, penalty, priorParameters, constrain, v=v,
        pass = 1, trace = trace, n=n, d=d, getTran=getTran, dependent=dependent) {
        
        g <- sample(1:(dim(mar$transformed)[[1]]), size = n, replace = TRUE)
        g <- mar$transformed[g, ]
        ok <- FALSE

        while (!ok) {
            for (j in 1:(dim(g)[[2]])){
				if (margins == "gumbel"){
					g[order(g[, j]), j] <- sort(-log(-log(runif(dim(g)[[1]]))))
				}
				else {
					u <- runif(nrow(g))
					g[order(g[, j]), j] <- sort(sign(u - .5) * log(1 - 2*abs(u - .5)))
				}
			}
            if (sum(g[, which] > dth) > 1){ ok <- TRUE }
        }
        
        g <- sapply(1:d, getTran, x = g, data = mar$data, margins=margins,
                    mod = mar$models, th = mar$mth, qu = mar$mqu)
                    
        dimnames(g)[[2]] <- names(mar$models)

        ggpd <- migpd(g, mth = mar$mth, 
					  penalty = penalty, priorParameters = priorParameters)

        gd <- mexDependence(ggpd, dth = dth, which = which, margins=margins, constrain=constrain, v=v)
        res <- list(GPD = coef(ggpd)[3:4, ],
                    dependence = gd$dependence$coefficients, 
                    Z = gd$dependence$Z,
                    Y = g)
                    
        if (pass == 1) {
            if (i%%trace == 0) {
                cat(paste(i, "replicates done\n"))
            }
        }
        res
    } # Close innerFun 

    res <- lapply(1:R, innerFun, x = x, which = which, dth = dth, margins=margins, 
        dqu = dqu, penalty = penalty, priorParameters = priorParameters, constrain=constrain, v=v,
        pass = 1, trace = trace, getTran=getTran, n=n, d=d, dependent=dependent)

    # Sometimes samples contain no extreme values. Need to have another pass or two
    if (nPass > 1) {
        for (pass in 2:nPass) {
            rerun <- sapply(res, function(x) any(sapply(x, function(x) any(is.na(x)))))
            wh <- !unlist(lapply(res, function(x) dim(x$Z)[[1]] > 0))
            rerun <- apply(cbind(rerun, wh), 1, any)
            if (sum(rerun) > 0) {
                cat("Pass", pass, ":", sum(rerun), "samples to rerun.\n")
                rerun <- (1:R)[rerun]
                res[rerun] <- lapply((1:R)[rerun], innerFun, 
                  x = x, which = which, dth = dth, dqu = dqu, margins=margins,
                  penalty = penalty, priorParameters = priorParameters, constrain=constrain, v=v,
                  pass = pass, trace = trace, getTran=getTran, n=n, d=d, dependent=dependent)
            }
        }
    }

    ans$boot <- res
    oldClass(ans) <- c("bootmex", "mex")
    ans
}

test.bootmex <- function(){ # this is a weak test - it tests the structure 
# of the output but not the correctness of the bootstrap coefficients; it will 
# also catch ERRORs (as opposed to FAILUREs) if the code breaks.  For strong 
# testing of this function, run test.predict.mex

  set.seed(20111007)
  smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
  wmarmod <- migpd(winter, mqu=.7,  penalty="none")

  mySdep <- mexDependence(smarmod, dqu=.7, which=1)
  myWdep <- mexDependence(wmarmod, dqu=.7, which=1)

  R <- 20

  mySboot <- bootmex(mySdep, R=R)
  myWboot <- bootmex(myWdep, R=R)

  checkEqualsNumeric(coef(mySdep)[[2]], mySboot$simpleDep, msg="bootmex: summer simpleDep from call with mar model")
  checkEqualsNumeric(coef(myWdep)[[2]], myWboot$simpleDep, msg="bootmex: winter simpleDep from call with mar model")

  checkEqualsNumeric(coef(smarmod), coef(mySboot$simpleMar), msg="bootmex: summer simpleMar from call with mar model")
  checkEqualsNumeric(coef(wmarmod), coef(myWboot$simpleMar), msg="bootmex: winter simpleMar from call with mar model")
  
  checkEqualsNumeric(dim(mySboot$boot[[1]]$Z)[2], dim(mySdep$dependence$Z)[2], msg="bootmex: summer dim of residuals from call with mar model")
  checkEqualsNumeric(dim(myWboot$boot[[1]]$Z)[2], dim(myWdep$dependence$Z)[2], msg="bootmex: winter dim of residuals from call with mar model")
  
  checkEqualsNumeric(dim(mySboot$boot[[1]]$dependence), dim(coef(mySdep)[[2]]), msg="bootmex: summer dim of coefficients from call with mar model")
  checkEqualsNumeric(dim(myWboot$boot[[1]]$dependence), dim(coef(myWdep)[[2]]), msg="bootmex: winter dim of coefficients from call with mar model")
  
  checkEqualsNumeric(R, length(mySboot$boot),msg="bootmex: number of bootstrap samples")
  checkEqualsNumeric(R, length(myWboot$boot),msg="bootmex: number of bootstrap samples")
  
  checkEqualsNumeric(dim(summer),dim(mySboot$boot[[1]]$Y),msg="bootmex: size of bootstrap data set")
  checkEqualsNumeric(dim(winter),dim(myWboot$boot[[5]]$Y),msg="bootmex: size of bootstrap data set")

  smexmod <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="gumbel",constrain=FALSE)
  wmexmod <- mex(winter, mqu=.7,  penalty="none", margins="gumbel",constrain=FALSE)

  mySboot <- bootmex(smexmod, R=R)
  myWboot <- bootmex(wmexmod, R=R)

  checkEqualsNumeric(coef(smexmod)[[2]], mySboot$simpleDep, msg="bootmex: summer simpleDep from call with mex model")
  checkEqualsNumeric(coef(smexmod)[[1]], coef(mySboot$simpleMar), msg="bootmex: summer simpleMar from call with mex model")

  checkEqualsNumeric(coef(wmexmod)[[2]], myWboot$simpleDep, msg="bootmex: winter simpleDep from call with mex model")
  checkEqualsNumeric(coef(wmexmod)[[1]], coef(myWboot$simpleMar), msg="bootmex: winter simpleMar from call with mex model")
  
  checkEqualsNumeric(R, length(mySboot$boot),msg="bootmex: number of bootstrap samples")
  checkEqualsNumeric(R, length(myWboot$boot),msg="bootmex: number of bootstrap samples")
  
  checkEqualsNumeric(dim(summer),dim(mySboot$boot[[1]]$Y),msg="bootmex: size of bootstrap data set")
  checkEqualsNumeric(dim(winter),dim(myWboot$boot[[5]]$Y),msg="bootmex: size of bootstrap data set")
  
  smexmod.1 <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="laplace",constrain=FALSE)
  smexmod.2 <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="laplace",constrain=TRUE,v=2)
  mySboot.1 <- bootmex(smexmod.1,R=R)
  mySboot.2 <- bootmex(smexmod.2,R=R)

  checkEqualsNumeric(coef(smexmod.1)[[2]], mySboot.1$simpleDep, msg="bootmex: summer simpleDep from call with mex model, constrain=FASLE")
  checkEqualsNumeric(coef(smexmod.1)[[1]], coef(mySboot.1$simpleMar), msg="bootmex: summer simpleMar from call with mex model, constrain=FASLE")
  checkEqualsNumeric(coef(smexmod.2)[[2]], mySboot.2$simpleDep, msg="bootmex: summer simpleDep from call with mex model, v=2")
  checkEqualsNumeric(coef(smexmod.2)[[1]], coef(mySboot.2$simpleMar), msg="bootmex: summer simpleMar from call with mex model, v=2")
  
# check execution of for 2-d data

  wavesurge.fit <- migpd(wavesurge,mq=.7)
  wavesurge.mex <- mexDependence(wavesurge.fit, dqu=0.8,which=1)
  R <- 20
  
  wavesurge.boot <- bootmex(wavesurge.mex,R=R)

  checkEqualsNumeric(dim(wavesurge.boot$boot[[1]]$Z)[2],1,msg="bootmex: execution for 2-d data")
  checkEquals(dimnames(wavesurge.boot$boot[[1]]$Z)[[2]],names(wavesurge)[2],msg="bootmex: execution for 2-d data")
  checkEqualsNumeric(length(wavesurge.boot$boot),R,msg="bootmex: execution for 2-d data")
}
