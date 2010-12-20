bootmex <- 
    # Bootstrap inference for a conditional multivaratiate extremes model.
function (x, which, R = 100, dth, dqu, nPass = 3, trace = 10) {
    theCall <- match.call()

    getgum <- function(i, x, data, mod, th, qu) {
        x <- c(x[, i])
        param <- mod[[i]]$coefficients
        th <- th[i]
        qu <- qu[i]
        data <- c(data[, i])
        res <- revGumbel(exp(-exp(-x)), data = data, th = th, 
            qu = qu, sigma = exp(param[1]), xi = param[2])
        res
    }

    if (class(x) == "mex"){
        if (!missing(which)){
            warning("which given, but already applied to 'mex' object. Using 'mex' value")
        }
        which <- x[[2]]$which
        if (!missing(dqu)){
            warning("dqu given, but already applied to 'mex' object. Using 'mex' value")
        }
        dqu <- x[[2]]$dqu
        x <- x[[1]]
    }
    else if (class(x) != "migpd"){
        stop("object should have class migpd")
    }
    else if (missing(which)) {
        cat("Missing 'which'. Conditioning on", dimnames(x$gumbel)[[2]][1], "\n")
        which <- 1
    }

    if (missing(dqu) & missing(dth)) {
        dqu <- x$mqu[which]
    }

    if (missing(dth)) {
        dth <- quantile(x$gumbel[, which], prob = dqu)
    }

    penalty <- x$penalty
    priorParameters <- x$priorParameters

    n <- dim(x$gumbel)[[1]]
    d <- dim(x$gumbel)[[2]]

    if (is.character(which)) {
        which <- (1:d)[dimnames(x$data)[[2]] == which]
    }
    dependent <- (1:d)[-which]

    innerFun <- function(i, x, which, dth, dqu, penalty, priorParameters, 
        pass = 1, trace = trace, n=n, d=d, getgum=getgum, dependent=dependent) {
        
        g <- sample(1:(dim(x$gumbel)[[1]]), size = n, replace = TRUE)
        g <- x$gumbel[g, ]
        ok <- FALSE

        while (!ok) {
            for (j in 1:(dim(g)[[2]])) g[order(g[, j]), j] <- sort(-log(-log(runif(dim(g)[[1]]))))
            if (sum(g[, which] > dth) > 1){ ok <- TRUE }
        }
        
        g <- sapply(1:d, getgum, x = g, data = x$data,
                    mod = x$models, th = x$mth, qu = x$mqu)
                    
        dimnames(g)[[2]] <- names(x$models)

        ggpd <- migpd(g, mth = x$mth, penalty = penalty, priorParameters = priorParameters)

        gd <- mexDependence(ggpd, dth = dth, which = which)
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

    dqu <- rep(dqu, d) 

    res <- lapply(1:R, innerFun, x = x, which = which, dth = dth, 
        dqu = dqu, penalty = penalty, priorParameters = priorParameters, 
        pass = 1, trace = trace, getgum=getgum, n=n, d=d, dependent=dependent)

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
                  x = x, which = which, dth = dth, dqu = dqu, 
                  penalty = penalty, priorParameters = priorParameters, 
                  pass = pass, trace = trace, getgum=getgum, n=n, d=d, dependent=dependent)
            }
        }
    }

    # Construct the object to be returned.
    ans <- list()
    ans$boot <- res
    ans$call <- theCall
    ans$dqu <- dqu
    ans$which <- which
    ans$R <- R
    ans$simpleMar <- x
    ans$simpleDep <- mexDependence(x, dth = dth, which)$coefficients
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

  checkEqualsNumeric(mySdep$coefficients, mySboot$simpleDep, msg="bootmex: summer simpleDep")
  checkEqualsNumeric(coef(smarmod), coef(mySboot$simpleMar), msg="bootmex: summer simpleMar")

  checkEqualsNumeric(myWdep$coefficients, myWboot$simpleDep, msg="bootmex: winter simpleDep")
  checkEqualsNumeric(coef(wmarmod), coef(myWboot$simpleMar), msg="bootmex: winter simpleMar")
  
  checkEqualsNumeric(R, length(mySboot$boot),msg="bootmex: number of bootstrap samples")
  checkEqualsNumeric(R, length(myWboot$boot),msg="bootmex: number of bootstrap samples")
  
  checkEqualsNumeric(dim(summer),dim(mySboot$boot[[1]]$Y),msg="cmxvBoot: size of bootstrap data set")
  checkEqualsNumeric(dim(winter),dim(myWboot$boot[[5]]$Y),msg="cmxvBoot: size of bootstrap data set")

  smarmod <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7)
  wmarmod <- mex(winter, mqu=.7,  penalty="none")

  mySboot <- bootmex(smarmod, R=R)
  myWboot <- bootmex(wmarmod, R=R)

  checkEqualsNumeric(coef(smarmod)[[2]], mySboot$simpleDep, msg="bootmex: summer simpleDep")
  checkEqualsNumeric(coef(smarmod)[[1]], coef(mySboot$simpleMar), msg="bootmex: summer simpleMar")

  checkEqualsNumeric(coef(wmarmod)[[2]], myWboot$simpleDep, msg="bootmex: winter simpleDep")
  checkEqualsNumeric(coef(wmarmod)[[1]], coef(myWboot$simpleMar), msg="bootmex: winter simpleMar")
  
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
