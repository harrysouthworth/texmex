gpd <- function(y, data, ...){
  if (!missing(data)) {
      y <- ifelse(deparse(substitute(y))== "substitute(y)", deparse(y),deparse(substitute(y)))
      y <- formula(paste(y, "~ 1"))
      y <- model.response(model.frame(y, data=data))
  } 
  UseMethod("gpd",y)
}

gpd.default <-
function (y, data, th, qu, phi = ~1, xi = ~1,
          penalty = "gaussian", prior = "gaussian",
          method = "optimize", cov="observed",
          start = NULL, priorParameters = NULL,
          maxit = 10000, trace=NULL,
          iter = 40500, burn=500, thin = 4,
          proposal.dist = c("gaussian", "cauchy"),
          jump.cov, jump.const, verbose=TRUE,...) {

    theCall <- match.call()

    ##################### Sort out penalties, priors, methods...

    if (!missing(penalty) & !missing(prior)){
        stop("specify one or neither of penalty and prior, not both")
    }
    if (!missing(penalty)){
        prior <- penalty
    }

    method <- casefold(method)
    prior <- casefold(prior)

    if (method %in% c("o", "opt", "optim", "optimize", "optimise")){
        method <- "o"
    }
    else if (method %in% c("s", "sim", "simulate")){
        method <- "s"
    }
    else {
        stop("method should be either 'optimize' or 'simulate'")
    }

    if (method == "s" & prior != "gaussian"){
        stop("only Gaussian prior can be used when simulating from the posterior")
    }

    ##################### Sort out trace
    if (method == "o"){
      if (!missing(trace)){
          otrace <- trace
      } else {
          otrace <- 0
      }
    } else{
      otrace <- 0
      if (missing(trace)){
         trace <- 1000
      }
    }
              
    ############################## Construct data to use...

    if (!missing(data)) {
      y <- deparse(substitute(y))
      y <- formula(paste(y, "~ 1"))
      y <- model.response(model.frame(y, data=data))
     
      if (!is.R() & length(as.character(phi)) == 2 & as.character(phi)[2] == "1"){
        X.phi <- matrix(rep(1, nrow(data)), ncol=1)
      } else {
        X.phi <- model.matrix(phi, data)
      }
    
      if (!is.R() & length(as.character(xi)) == 2 & as.character(xi)[2] == "1"){
        X.xi <- matrix(rep(1, nrow(data)), ncol=1)
      } else {
        X.xi <- model.matrix(xi, data)
      }
    } else {
      if (length(as.character(phi)) == 2 & as.character(phi)[2] == "1"){  
        X.phi <- matrix(ncol = 1, rep(1, length(y)))
      } else {
        X.phi <- model.matrix(phi)
      }
      if (length(as.character(xi)) == 2 & as.character(xi)[2] == "1"){
        X.xi <- matrix(ncol = 1, rep(1, length(y)))
      } else {
        X.xi <- model.matrix(xi)
      }
    }
    if (missing(th)) {
        th <- quantile(y, qu)
    }
    X.phi <- X.phi[y > th, ]
    X.xi <- X.xi[y > th, ]    
    rate <- mean(y > th)

    allY <- y
    y <- y[y > th]
    if (!is.matrix(X.phi)) {
        X.phi <- matrix(X.phi, ncol = 1)
    }
    if (!is.matrix(X.xi)) {
        X.xi <- matrix(X.xi, ncol = 1)
    }
    if (length(y) == 0) {
        stop("No observations over threshold")
    }

    ###################### Check and sort out prior parameters...

    if (prior %in% c("quadratic", "gaussian")) {
        if (is.null(priorParameters)) {
            priorParameters <- list(rep(0, ncol(X.phi) + ncol(X.xi)), 
                diag(rep(10^4, ncol(X.phi) + ncol(X.xi))))
        }
        if (length(priorParameters) != 2 | !is.list(priorParameters)) {
            stop("For Gaussian prior or quadratic penalty, priorParameters should be a list of length 2, the second element of which should be a symmetric (covariance) matrix")
        }
    } else if (prior %in% c("lasso", "l1", "laplace")) {
        if (is.null(priorParameters)) {
            priorParameters <- list(rep(0, ncol(X.phi) + ncol(X.xi)), 
                diag(rep(10^(-4), ncol(X.phi) + ncol(X.xi))))
        }
        if (length(priorParameters) != 2 | !is.list(priorParameters)) {
            stop("For Laplace prior or L1 or Lasso penalty, priorParameters should be a list of length 2, the second element of which should be a diagonal (precision) matrix")
        }
        if (!is.matrix(priorParameters[[2]])) {
            priorParameters[[2]] <- diag(rep(priorParameters[[2]], 
                ncol(X.phi) + ncol(X.xi)))
        }
        if (!all(priorParameters[[2]] == diag(diag(priorParameters[[2]])))) {
            warning("some off-diagonal elements of the covariance are non-zero. Only the diagonal is used in penalization")
        }
    }

    #### If priorParameters given but of wrong dimension, kill
    if (!is.null(priorParameters)) {
        dimp <- ncol(X.phi) + ncol(X.xi)
        if (length(priorParameters[[1]]) != dimp) {
            stop("wrong number of parameters in prior (doesn't match phi and xi formulas)")
        }
        else if (length(diag(priorParameters[[2]])) != dimp) {
            stop("wrong dimension of prior covariance (doesn't match phi and xi formulas)")
        }
    }

    ################################## Do the optimization....

    o <- gpdFit(y=y, th=th, X.phi=X.phi, X.xi=X.xi, penalty=prior, start=start, 
                 hessian = cov == "numeric",
                 priorParameters = priorParameters, maxit = maxit, trace = otrace)


    if (o$convergence != 0) {
        warning("Non-convergence in gpd")
    }

    ################################## If method = "optimize", construct object and return...

    if( cov == "observed" ){
      o$hessian <- NULL
    }
    o$threshold <- th
    o$penalty <- prior
    o$coefficients <- o$par
    names(o$coefficients) <- c(paste("phi:", colnames(X.phi)), 
                               paste("xi:", colnames(X.xi)))

    o$formulae <- list(phi=phi, xi=xi)
    o$par <- NULL
    o$rate <- rate
    o$call <- theCall
    o$y <- y
    o$X.phi <- X.phi
    o$X.xi <- X.xi
	  o$priorParameters <- priorParameters
	  if (missing(data)) {
        data <- allY
    }
    o$data <- data
                     
    fittedScale <- fittedGPDscale(o)
    fittedShape <- fittedGPDshape(o)
    scaledY <- fittedShape * (o$y - o$threshold) / fittedScale
    o$residuals <- 1/fittedShape * log(1 + scaledY) # Standard exponential
    if(any(fittedShape < -0.5)){
      warning("fitted shape parameter xi < -0.5\n")
    }
    o$loglik <- -o$value
    o$value <- NULL
    o$counts <- NULL
    oldClass(o) <- "gpd"

    if (cov == "numeric") {
      o$cov <- solve(o$hessian)
    } else if (cov == "observed") {
      o$cov <- solve(info.gpd(o))
    } else {
      stop("cov must be either 'numeric' or 'observed'")
    }

    o$se <- sqrt(diag(o$cov))
    if (method == "o"){
		  o
    }

    ################################# Simulate from posteriors....

    else { # Method is "simulate"...
      proposal.fn <- switch(match.arg(proposal.dist),
                            gaussian=rmvnorm,
                            cauchy=.rmvcauchy,
                            function () {stop("Bad proposal distribution")})
      if (missing(jump.const)){
        jump.const <- (2.4/sqrt(ncol(X.phi) + ncol(X.xi)))^2
      }
      u <- th

      prior <- .make.mvn.prior(priorParameters)

############################# Define log-likelihood

      gpd.log.lik <- .make.gpd.loglikelihood(y, u, X.phi, X.xi)

      log.lik <- function(param) {
        gpd.log.lik(param) + prior(param)
      }

      # Need to check for convergence failure here. Otherwise, end up simulating
      # proposals from distribution with zero variance in 1 dimension.
       checkNA <- any(is.na(sqrt(diag(o$cov))))
      if (checkNA){
        stop("MAP estimates have not converged. Cannot proceed. Try a different prior" )
      }

      res <- matrix(ncol=ncol(X.phi) + ncol(X.xi), nrow=iter)
      res[1,] <- if (missing(start)) { o$coefficients } else { start }


      if (!exists(".Random.seed")){ runif(1)  }
      seed <- .Random.seed # Retain and add to output

      cov <- if (missing(jump.cov)) { o$cov } else { jump.cov }
                                        # create proposals en bloc
      proposals <- proposal.fn(iter,
                               double(length(o$coefficients)),
                               cov*jump.const)
      last.cost <- log.lik(res[1,])
      if (!is.finite(last.cost)) {
        stop("Start is infeasible.")
      }
######################## Run the Metropolis algorithm...
      acc <- 0
      for(i in 2:iter){
        if( verbose){
          if( i %% trace == 0 ) cat(i, " steps taken\n" )
        }
        prop <- proposals[i - 1,] + res[i - 1,]
        top <- log.lik(prop)
        delta <- top - last.cost
        if (is.finite(top) && ((delta >= 0) ||
                               (runif(1) <= exp(delta)))) {
          res[i, ] <- prop
          last.cost <- top
          acc <- 1 + acc
        } else {
          res[ i , ] <- res[i-1,]
        }
      } # Close for(i in 2:iter

      acc <- acc / iter
      if (acc < .1) {
        warning("Acceptance rate in Metropolis algorithm is low.")
      }

      if (trace < iter) {
        if(verbose) {
          cat("Acceptance rate:", round( acc , 3 ) , "\n")
        }
      }

      if (missing(data)) {
        data <- allY
      }
      res <- list(call=theCall, threshold=u , map = o,
                  burn = burn, thin = thin,
                  chains=res, y=y, data=data,
                  X.phi = X.phi, X.xi = X.xi,
                  acceptance=acc, seed=seed)

      oldClass(res) <- "bgpd"
      res <- thinAndBurn(res)
      res
    } # Close else

}

test.gpd <- function(){
  tol <- 0.01

###################################################################
# 1.3 Reproduce loglik, parameter estimates and covariance on page 85
#    of Coles. Will not be exact because fitting routines differ:
#    gpd works with phi=log(sigma). Can only reproduce cell [2,2]
#    of covariance.

  cparas <- c(7.44, 0.184)
  cse <- c(0.958432, 0.101151)

  ccov <- matrix(c(.9188, -.0655, -.0655, .0102), nrow=2)
  cloglik <- -485.1

  mod <- gpd(rain, th=30, penalty="none")
  mod.coef <- coef(mod)

  mod.coef[1] <- exp(mod.coef[1])
  names(mod.coef)[1] <- "sigma"

  mod.loglik <- mod$loglik
  mod.cov22 <- mod$cov[2, 2]

  checkEqualsNumeric(cparas, mod.coef, tolerance = tol,
                     msg="gpd: parameter ests page 85 Coles")
  checkEqualsNumeric(cse[2], mod$se[2], tolerance = tol,
                     msg="gpd: standard errors page 85 Coles")
  checkEqualsNumeric(ccov[2,2], mod$cov[2, 2], tolerance = tol,
                     msg="gpd: Covariance page 85 Coles")
  checkEqualsNumeric(cloglik, mod.loglik, tolerance = tol,
                     msg="gpd: loglik page 85 Coles")

###################################################################
#   Logical checks on the effect of Gaussian penalization. The smaller the
#    variance, the more the parameter should be drawn towards the
#    mean.

# 2.1 Tests for xi being drawn to 0

  gp1 <- list(c(0, 0), diag(c(10^4, .25)))
  gp2 <- list(c(0, 0), diag(c(10^4, .05)))

  mod1 <- gpd(rain, th=30, priorParameters=gp1)
  mod2 <- gpd(rain, th=30, priorParameters=gp2)

  checkTrue(coef(mod)[2] > coef(mod1)[2],
            msg="gpd: Gaussian penalization xi being drawn to 0")
  checkTrue(coef(mod1)[2] > coef(mod2)[2],
            msg="gpd: Gaussian penalization xi being drawn to 0")

# 2.2 Tests for phi being drawn to 0

  gp3 <- list(c(0, 0), diag(c(1, 10^4)))
  gp4 <- list(c(0, 0), diag(c(.1, 10^4)))

  mod3 <- gpd(rain, th=30, priorParameters=gp3)
  mod4 <- gpd(rain, th=30, priorParameters=gp4)

  checkTrue(coef(mod)[1] > coef(mod3)[1],
            msg="gpd: Gaussian penalization phi being drawn to 0")
  checkTrue(coef(mod3)[1] > coef(mod4)[1],
            msg="gpd: Gaussian penalization phi being drawn to 0")
  
# 2.3 Tests for xi being drawn to 1
  gp5 <- list(c(0, 1), diag(c(10^4, .25)))
  gp6 <- list(c(0, 1), diag(c(10^4, .05)))

  mod5 <- gpd(rain, th=30, priorParameters=gp5)
  mod6 <- gpd(rain, th=30, priorParameters=gp6)

  checkTrue(1 - coef(mod)[2] > 1 - coef(mod5)[2],
            msg="gpd: Gaussian penalization xi being drawn to 1")
  checkTrue(1 - coef(mod1)[2] > 1 - coef(mod6)[2],
            msg="gpd: Gaussian penalization xi being drawn to 1")
  
# 2.4 Tests for phi being drawn to 4 (greater than mle for phi)

  gp7 <- list(c(4, 0), diag(c(1, 10^4)))
  gp8 <- list(c(4, 0), diag(c(.1, 10^4)))

  mod7 <- gpd(rain, th=30, priorParameters=gp7)
  mod8 <- gpd(rain, th=30, priorParameters=gp8)

  checkTrue(4 - coef(mod)[1] > 4 - coef(mod7)[1],
            msg="gpd: Gaussian penalization phi being drawn to 4")
  checkTrue(4 - coef(mod3)[1] > 4 - coef(mod8)[1],
            msg="gpd: Gaussian penalization phi being drawn to 4")
  
###################################################################
#   Logical checks on the effect of penalization using lasso or L1 penalization. The smaller the
#    variance, the more the parameter should be drawn towards the
#    mean.

# 2a.1 Tests for xi being drawn to 0

  gp1 <- list(c(0, 0), solve(diag(c(10^4, .25))))
  gp2 <- list(c(0, 0), solve(diag(c(10^4, .05))))

  mod1 <- gpd(rain, th=30, priorParameters=gp1, penalty="lasso")
  mod2 <- gpd(rain, th=30, priorParameters=gp2, penalty="lasso")

  checkTrue(coef(mod)[2] > coef(mod1)[2],
            msg="gpd: lasso penalization xi being drawn to 0")
  checkTrue(coef(mod1)[2] > coef(mod2)[2],
            msg="gpd: lasso penalization xi being drawn to 0")

# 2a.2 Tests for phi being drawn to 0

  gp3 <- list(c(0, 0), solve(diag(c(1, 10^4))))
  gp4 <- list(c(0, 0), solve(diag(c(.1, 10^4))))

  mod3 <- gpd(rain, th=30, priorParameters=gp3, penalty="lasso")
  mod4 <- gpd(rain, th=30, priorParameters=gp4, penalty="lasso")

  checkTrue(coef(mod)[1] > coef(mod3)[1],
            msg="gpd: lasso penalization phi being drawn to 0")
  checkTrue(coef(mod3)[1] > coef(mod4)[1],
            msg="gpd: lasso penalization phi being drawn to 0")
  
# 2a.3 Tests for xi being drawn to 1
  gp5 <- list(c(0, 1), solve(diag(c(10^4, .25))))
  gp6 <- list(c(0, 1), solve(diag(c(10^4, .05))))

  mod5 <- gpd(rain, th=30, priorParameters=gp5, penalty="lasso")
  mod6 <- gpd(rain, th=30, priorParameters=gp6, penalty="lasso")

  checkTrue(1 - coef(mod)[2] > 1 - coef(mod5)[2],
            msg="gpd: lasso penalization xi being drawn to 1")
  checkTrue(1 - coef(mod1)[2] > 1 - coef(mod6)[2],
            msg="gpd: lasso penalization xi being drawn to 1")
  
# 2a.4 Tests for phi being drawn to 4 (greater than mle for phi)

  gp7 <- list(c(4, 0), solve(diag(c(1, 10^4))))
  gp8 <- list(c(4, 0), solve(diag(c(.1, 10^4))))

  mod7 <- gpd(rain, th=30, priorParameters=gp7, penalty="lasso")
  mod8 <- gpd(rain, th=30, priorParameters=gp8, penalty="lasso")

  checkTrue(4 - coef(mod)[1] > 4 - coef(mod7)[1],
            msg="gpd: lasso penalization phi being drawn to 4")
  checkTrue(4 - coef(mod3)[1] > 4 - coef(mod8)[1],
            msg="gpd: lasso penalization phi being drawn to 4")

########################################################
# Tests on including covariates. Once more, gpd.fit in ismev
# works with sigma inside the optimizer, so we need to tolerate
# some differences and standard errors might be a bit out.

# 3.0 Reproduce Coles, page 119. Reported log-likelihood is -484.6.

  rtime <- (1:length(rain))/1000
  d <- data.frame(rainfall = rain, time=rtime)

  mod <- gpd(rainfall, th=30, data=d, phi= ~ time, penalty="none")

  checkEqualsNumeric(-484.6, mod$loglik, tolerance = tol,
                     msg="gpd: loglik Coles page 119")
  
####################################################################
# 3.1 Use liver data, compare with ismev. 
#     These are not necessarily sensible models!
#     Start with phi alone.

  mod <- gpd(ALT.M, qu=.7, data=liver,
           phi = ~ ALT.B + dose, xi = ~1,
           penalty="none", cov="observed")

  m <- model.matrix(~ ALT.B + dose, liver)

  ismod <- texmex:::.ismev.gpd.fit(liver$ALT.M,
                                   threshold=quantile(liver$ALT.M, .7),
                                   ydat = m, sigl=2:ncol(m),
                                   siglink=exp, show=FALSE)

  checkEqualsNumeric(ismod$mle,
                     coef(mod),
                     tolerance=tol,
                     msg="gpd: covariates in phi only, point ests")

# SEs for phi will not be same as for sigma, but we can test xi
  checkEqualsNumeric(ismod$se[length(ismod$se)],
                     mod$se[length(mod$se)],
                     tolerance=tol,
                     msg="gpd: covariates in phi only, standard errors")

######################################################################
# 3.2 Test xi alone.
  mod <- gpd(log(ALT.M / ALT.B), qu=.7, data=liver,
           phi = ~1, xi = ~ ALT.B + dose,
           penalty="none")

  m <- model.matrix(~ ALT.B + dose, liver)

  ismod <- texmex:::.ismev.gpd.fit(log(liver$ALT.M / liver$ALT.B), 
                                   threshold=quantile(log(liver$ALT.M / liver$ALT.B), .7), 
                                   ydat = m, shl=2:ncol(m), show=FALSE)
  mco <- coef(mod)
  mco[1] <- exp(mco[1])

  checkEqualsNumeric(ismod$mle,
                     mco,
                     tolerance=tol,
                     msg="gpd: covariates in xi only: point ests")
# SEs for phi will not be same as for sigma, but we can test xi
  checkEqualsNumeric(ismod$se[-1],
                     mod$se[-1],
                     tolerance = tol,
                     msg="gpd: covariates in xi only: standard errors")

######################################################################
# 3.3 Test phi & xi simultaneously. Use simulated data.

  set.seed(25111970)
  
  makeData <- function(a,b,n=500,u=10)
  # lengths of a and b should divide n exactly
  # returns data set size 2n made up of uniform variates (size n) below threshold u and 
  # gpd (size n) with scale parameter exp(a) and shape b above threshold u
  {
    gpd <- rgpd(n,exp(a),b,u=u)
    unif <- runif(n,u-10,u)
    as.data.frame(cbind(a=a,b=b,y=c(gpd,unif)))
  }

  mya <- seq(0.1,1,len=10)
  myb <- rep(c(-0.2,0.2),each=5)
  data <- makeData(mya,myb)
  m <- model.matrix(~ a+b, data)
  
  mod <- gpd(y,qu=0.7,data=data,phi=~a,xi=~b,penalty="none")
  ismod <- texmex:::.ismev.gpd.fit(data$y,
                                   threshold=quantile(data$y,0.7),
                                   ydat=m,shl=3,sigl=2,
                                   siglink=exp,
                                   show=FALSE)

  checkEqualsNumeric(ismod$mle,
                     coef(mod),
                     tolerance = tol,
                     msg="gpd: covariates in phi and xi: point ests")
  checkEqualsNumeric(ismod$se,
                     sqrt(diag(mod$cov)),
                     tolerance = tol,
                     msg="gpd: covariates in phi and xi: std errs")

####################################################################
# Check that using priors gives expected behaviour when covariates are included.

# 2.1 Tests for xi being drawn to 0

  myb <- rep(c(0.5,1.5),each=5)
  data <- makeData(a=1,b=myb,n=3000)
  
  gp1 <- list(c(0, 0, 0), diag(c(10^4, 0.25, 0.25)))
  gp2 <- list(c(0, 0, 0), diag(c(10^4, 0.25, 0.01)))

  mod0 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,penalty="none")
  mod1 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp1)
  mod2 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp2)

  checkTrue(all(abs(coef(mod0)[2:3]) > abs(coef(mod1)[2:3])),
            msg="gpd: with covariates, xi drawn to zero")
  checkTrue(abs(coef(mod1)[3]) > abs(coef(mod2)[3]),
            msg="gpd: with covariates, xi drawn to zero")

# 2.2 Tests for phi being drawn to 0

  # HS. Changed a to mya due to scoping problems in S+. The issue is very general
  # and affects (for example) lm(~a, data, method="model.frame"), so it's kind of
  # by design.
  mya <- seq(0.1,1,len=10)
  data <- makeData(-3 + mya,b=-0.1,n=3000)
  data$a <- rep(mya, len=nrow(data))

  gp4 <- list(c(0, 0, 0), diag(c(1, 1, 10^4)))
  gp5 <- list(c(0, 0, 0), diag(c(0.1, 0.1, 10^4)))

  mod3 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,penalty="none")
  mod4 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp4)
  mod5 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp5)

  checkTrue(all(abs(coef(mod3)[1:2]) > abs(coef(mod4)[1:2])),
            msg="gpd: with covariates, phi being drawn to 0")
  checkTrue(all(abs(coef(mod4)[1:2]) > abs(coef(mod5)[1:2])),
            msg="gpd: with covariates, phi being drawn to 0")

# 2.3 Tests for xi being drawn to 2
  myb <- rep(c(-0.5,0.5),each=5)
  data <- makeData(a=1,b=myb,n=3000)
 
  gp7 <- list(c(0, 2, 2), diag(c(10^4, 0.25, 0.25)))
  gp8 <- list(c(0, 2, 2), diag(c(10^4, 0.05, 0.05)))

  mod6 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,penalty="none")
  mod7 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp7)
  mod8 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp8)

  checkTrue(all(abs(2 - coef(mod6)[2:3]) > abs(2 - coef(mod7)[2:3])),
            msg="gpd: with covariates, xi drawn to 2")
  checkTrue(all(abs(2 - coef(mod7)[2:3]) > abs(2 - coef(mod8)[2:3])),
            msg="gpd: with covariates, xi drawn to 2")

# 2.4 Tests for phi being drawn to 4 

  mya <- seq(0.1,1,len=10)
  data <- makeData(2 + mya,b=-0.1,n=3000)
  data$a <- rep(mya, len=nrow(data))
  
  gp10 <- list(c(0, 4, 0), diag(c(10^4, 1,   10^4)))
  gp11 <- list(c(0, 4, 0), diag(c(10^4, 0.1, 10^4)))

  mod9 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,penalty="none")
  mod10 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp10)
  mod11 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp11)

  checkTrue(abs(4 - coef(mod9)[2])  > abs(4 - coef(mod10)[2]),
            msg="gpd: with covariates, phi drawn to 4")
  checkTrue(abs(4 - coef(mod10)[2]) > abs(4 - coef(mod11)[2]),
            msg="gpd: with covariates, phi drawn to 4")

#*************************************************************
  postSum <- function(x){
    t(apply(x$param, 2, function(o){ c(mean=mean(o), se=sd(o)) }))
  } 

#************************************************************* 
# 4.1. Test reproducibility
  set.seed(20101110)
  save.seed <- .Random.seed

  set.seed(save.seed)
  bmod <- gpd(ALT.M, data=liver,
              th=quantile(liver$ALT.M, .7),
              iter=1000, thin=1, verbose=FALSE, method="sim")
  
  set.seed(save.seed)
  bmod2 <- gpd(ALT.M, data=liver,
               th=quantile(liver$ALT.M, .7),
               iter=1000, thin=1, verbose=FALSE, method="sim")
  
  checkEqualsNumeric(bmod$param, bmod2$param,
                     msg="gpd: test simulation reproducibility 1")

  set.seed(bmod$seed)
  bmod3 <- gpd(ALT.M, data=liver,
               th=quantile(liver$ALT.M, .7),
               iter=1000, thin=1, verbose=FALSE, method="sim")
  checkEqualsNumeric(bmod$param, bmod3$param,
                     msg="gpd: test simulation reproducibility 2")

#*************************************************************  
# 4.2. Logical test of burn-in

  checkEqualsNumeric(nrow(bmod$chains) - bmod$burn, nrow(bmod$param),
                     msg="gpd: Logical test of burn-in 1")

  iter <- sample(500:1000,1)
  burn <- sample(50,1)
  bmod2 <- gpd(ALT.M, data=liver, th=quantile(liver$ALT.M, .7),
                iter=iter, burn=burn, thin=1, verbose=FALSE, method="sim")

  checkEqualsNumeric(iter-burn, nrow(bmod2$param),
                     msg="gpd: Logical test of burn-in 2")

#*************************************************************
# 4.3. Logical test of thinning

  thin <- 0.5
  iter <- 1000
  bmod <- gpd(ALT.M, data=liver, th=quantile(liver$ALT.M, .7),
               iter=iter, thin = thin,verbose=FALSE, method="sim")

  checkEqualsNumeric((nrow(bmod$chains) - bmod$burn) * thin, nrow(bmod$param),
                     msg="gpd: Logical test of thinning 1")

  thin <- 2
  iter <- 1000
  bmod <- gpd(ALT.M, data=liver, th=quantile(liver$ALT.M, .7),
               iter=iter, thin = thin, verbose=FALSE, method="sim")

  checkEqualsNumeric((nrow(bmod$chains) - bmod$burn) / thin, nrow(bmod$param),
                     msg="gpd: Logical test of thinning 1")

}

