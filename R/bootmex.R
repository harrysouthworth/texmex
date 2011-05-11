bootmex <- 
    # Bootstrap inference for a conditional multivaratiate extremes model.
function (x, R = 100, nPass = 3, trace = 10) {
    theCall <- match.call()

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

    if (class(x) != "mexDependence"){
      stop("object must be of type 'mexDependence'")
    }
    which <- x$which
    dqu <- x$dqu
    dth <- x$dth
    margins <- x$migpd$margins        
    penalty <- x$migpd$penalty
    priorParameters <- x$migpd$priorParameters

    n <- dim(x$migpd$transformed)[[1]]
    d <- dim(x$migpd$transformed)[[2]]
    dqu <- rep(dqu, d) 
    dependent <- (1:d)[-which]
    
    ans$simpleDep <- x$coefficients
    ans$dqu <- dqu
    ans$which <- which
    ans$R <- R
    ans$simpleMar <- x$migpd
    ans$margins <- margins

    innerFun <- function(i, x, which, dth, dqu, margins, penalty, priorParameters, 
        pass = 1, trace = trace, n=n, d=d, getTran=getTran, dependent=dependent) {
        
        g <- sample(1:(dim(x$migpd$transformed)[[1]]), size = n, replace = TRUE)
        g <- x$migpd$transformed[g, ]
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
        
        g <- sapply(1:d, getTran, x = g, data = x$migpd$data, margins=margins,
                    mod = x$migpd$models, th = x$migpd$mth, qu = x$migpd$mqu)
                    
        dimnames(g)[[2]] <- names(x$migpd$models)

        ggpd <- migpd(g, mth = x$migpd$mth, 
					  penalty = penalty, priorParameters = priorParameters)

        gd <- mexDependence(ggpd, dth = dth, which = which, margins=margins, start = ans$simpleDep[c(1:2,5:6),])
        res <- list(GPD = coef(ggpd)[3:4, ],
                    dependence = gd$coefficients, 
                    Z = gd$Z,
                    Y = g)
                    
        if (pass == 1) {
            if (i%%trace == 0) {
                cat(paste(i, "replicates done\n"))
            }
        }
        res
    } # Close innerFun 



    res <- lapply(1:R, innerFun, x = x, which = which, dth = dth, margins=margins, 
        dqu = dqu, penalty = penalty, priorParameters = priorParameters, 
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
                  penalty = penalty, priorParameters = priorParameters, 
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

  smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
  wmarmod <- migpd(winter, mqu=.7,  penalty="none")

  mySdep <- mexDependence(smarmod, dqu=.7)
  myWdep <- mexDependence(wmarmod)

  R <- 20
 
  mySboot <- bootmex(smarmod, R=R, which=1, dqu=.7)
  myWboot <- bootmex(wmarmod, R=R, which=1, dqu=.7)

  checkEqualsNumeric(mySdep$coefficients, mySboot$simpleDep, msg="bootmex: summer simpleDep from call with mar model")
  checkEqualsNumeric(coef(smarmod), coef(mySboot$simpleMar), msg="bootmex: summer simpleMar from call with mar model")

  checkEqualsNumeric(myWdep$coefficients, myWboot$simpleDep, msg="bootmex: winter simpleDep from call with mar model")
  checkEqualsNumeric(coef(wmarmod), coef(myWboot$simpleMar), msg="bootmex: winter simpleMar from call with mar model")
  
  checkEqualsNumeric(R, length(mySboot$boot),msg="bootmex: number of bootstrap samples")
  checkEqualsNumeric(R, length(myWboot$boot),msg="bootmex: number of bootstrap samples")
  
  checkEqualsNumeric(dim(summer),dim(mySboot$boot[[1]]$Y),msg="cmxvBoot: size of bootstrap data set")
  checkEqualsNumeric(dim(winter),dim(myWboot$boot[[5]]$Y),msg="cmxvBoot: size of bootstrap data set")

  smexmod <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="gumbel")
  wmexmod <- mex(winter, mqu=.7,  penalty="none", margins="gumbel")

  mySboot <- bootmex(smexmod, R=R)
  myWboot <- bootmex(wmexmod, R=R)

  checkEqualsNumeric(coef(smexmod)[[2]], mySboot$simpleDep, msg="bootmex: summer simpleDep from call with mex model")
  checkEqualsNumeric(coef(smexmod)[[1]], coef(mySboot$simpleMar), msg="bootmex: summer simpleMar from call with mex model")

  checkEqualsNumeric(coef(wmexmod)[[2]], myWboot$simpleDep, msg="bootmex: winter simpleDep from call with mex model")
  checkEqualsNumeric(coef(wmexmod)[[1]], coef(myWboot$simpleMar), msg="bootmex: winter simpleMar from call with mex model")
  
  checkEqualsNumeric(R, length(mySboot$boot),msg="bootmex: number of bootstrap samples")
  checkEqualsNumeric(R, length(myWboot$boot),msg="bootmex: number of bootstrap samples")
  
  checkEqualsNumeric(dim(summer),dim(mySboot$boot[[1]]$Y),msg="cmxvBoot: size of bootstrap data set")
  checkEqualsNumeric(dim(winter),dim(myWboot$boot[[5]]$Y),msg="cmxvBoot: size of bootstrap data set")

# check execution of for 2-d data

  wavesurge.fit <- migpd(wavesurge,mq=.7)
  R <- 20
  wavesurge.boot <- bootmex(wavesurge.fit,which=1,R=R,dqu=0.8)

  checkEqualsNumeric(dim(wavesurge.boot$boot[[1]]$Z)[2],1,msg="bootmex: execution for 2-d data")
  checkEquals(dimnames(wavesurge.boot$boot[[1]]$Z)[[2]],names(wavesurge)[2],msg="bootmex: execution for 2-d data")
  checkEqualsNumeric(length(wavesurge.boot$boot),R,msg="bootmex: execution for 2-d data")
  
}
