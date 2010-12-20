gpd <- 
function (y, data, th, qu, phi = ~1, xi = ~1,
          penalty = "gaussian", prior = "gaussian", method = "optimize",
          start = NULL, priorParameters = NULL, maxit = 10000, trace=NULL,
          iter = 10500, burn=500, thin = 1, jump.const, verbose=TRUE) {

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
        }
        else {
            otrace <- 0
        }
    }
    else{
       otrace <- 0
       if (missing(trace)){
           trace <- 1000
       }
    }
              
    ############################## Construct data to use...

    if (!missing(data)) {
      y <- deparse(substitute(y))
    }
    if (!missing(data)) {
        y <- formula(paste(y, "~ 1"))
        y <- model.response(model.frame(y, data=data))
        if (missing(th)) {
            th <- quantile(y, qu)
        }
        X.phi <- model.matrix(phi, data)
        X.xi <- model.matrix(xi, data)
        X.phi <- X.phi[y > th, ]
        X.xi <- X.xi[y > th, ]
    }
    else {
        if (missing(th)) {
            th <- quantile(y, qu)
        }
        X.phi <- matrix(ncol = 1, rep(1, length(y)))
        X.xi <- matrix(ncol = 1, rep(1, length(y)))
        X.phi <- X.phi[y > th, ]
        X.xi <- X.xi[y > th, ]
    }
    rate <- mean(y > th)


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
    }
    else if (prior %in% c("lasso", "l1", "laplace")) {
        if (is.null(priorParameters)) {
            priorParameters <- list(rep(0, ncol(X.phi) + ncol(X.xi)), 
                diag(rep(10^(-4), ncol(X.phi) + ncol(X.xi))))
        }
        if (length(priorParameters) != 2 | !is.list(priorParameters)) {
            stop("For Laplace prior, priorParameters should be a list of length 2, the second element of which should be a diagonal matrix")
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

    o <- gpd.fit(y=y, th=th, X.phi=X.phi, X.xi=X.xi, penalty=prior,
                 priorParameters = priorParameters, maxit = maxit, trace = otrace)


    if (o$convergence != 0) {
        warning("Non-convergence in gpd")
    }

    ################################## If method = "optimize", construct object and return...

    if (method == "o"){
        o$cov <- solve(o$hessian)
        o$hessian <- NULL
        o$se <- sqrt(diag(o$cov))
        o$threshold <- th
        o$penalty <- prior
        o$coefficients <- o$par
        names(o$coefficients) <- c(paste("phi:", colnames(X.phi)), 
                                   paste("xi:", colnames(X.xi)))
        o$par <- NULL
        o$rate <- rate
        o$call <- theCall
        o$y <- y
        o$X.phi <- X.phi
        o$X.xi <- X.xi
        if (missing(data)) {
            data <- NULL
        }
        o$data <- data

        o$loglik <- -o$value
        o$value <- NULL
        o$counts <- NULL
        oldClass(o) <- "gpd"
        o
    } # Close if (method == "o"

    ################################# Simulate from posteriors....

    else { # Method is "simulate"...
   
        if (missing(jump.const)){
            jump.const <- (2.4/sqrt(ncol(X.phi) + ncol(X.xi)))^2
        }

        u <- th
	
        # Next line seems redundant, but leave in case want to allow other priors in future
        prior <- dmvnorm

        gpdlik <- # Positive loglikelihood
        function(param, data, X.phi, X.xi){
            phi <- colSums(param[1:ncol(X.phi)] * t(X.phi))
            xi <- colSums(param[(ncol(X.phi) + 1):(ncol(X.phi) + ncol(X.xi))] * t(X.xi))

            n <- length(data)
            data <- xi * data / exp(phi) + 1
            if(any(data <= 0)){ #| ( xi <= -.50 ) ) -(10^10)
                -(10^10)
            }
            else  {
                - sum(phi) - sum(log(data) * (1/xi + 1))
            }
        } # Close gpdlik <- function



        # Need to check for convergence failure here. Otherwise, end up simulating
        # proposals from distribution with zero variance in 1 dimension.
        checkNA <- any(is.na(sqrt(diag(o$cov))))
        if ( checkNA ){
            stop( "MAP estimates have not converged. Cannot proceed. Try a different prior" )
        }

        res <- matrix( ncol=ncol(X.phi) + ncol(X.xi), nrow=iter )

        if (missing(start)) {
            res[1, ] <- o$par
        }
	    else { 
            res[1, ] <- start
        }

        if (!exists(".Random.seed")){ runif(1)  }
        seed <- .Random.seed # Retain and add to output


        cov <- solve(o$hessian)
        ######################## Run the Metropolis algorithm...
        for(i in 2:iter){
            if( verbose){
                if( i %% trace == 0 ) cat(i, " steps taken\n" )
            }
        
            if ( exists("is.R") && is.function(is.R) && is.R() ){
                prop <- c(rmvnorm( 1 , mean = o$par, sigma = cov * jump.const ))
            }
            else {
                prop <- c(rmvnorm( 1 , mean = o$par, cov = cov * jump.const ))
            }
            bottom <- prior(res[i - 1,], priorParameters[[1]], priorParameters[[2]], log=TRUE ) +
                      gpdlik(res[ i - 1 , ] , y-u, X.phi, X.xi)
            top <- prior(prop, priorParameters[[1]], priorParameters[[2]], log=TRUE) +
                   gpdlik(prop, y-u, X.phi, X.xi)
            r <- exp(top - bottom)
            res[ i , ] <- if (runif(1) < r ){ prop }
                          else{ res[i -1, ] }
        } # Close for(i in 2:iter

        acc <- length(unique(res[,1])) / iter
        if (acc < .1){ warning("Acceptance rate in Metropolis algorithm is low.") }
        if (trace < iter){
            if(verbose) {
                cat("Acceptance rate:", round( acc , 3 ) , "\n")
            }
        }

        res <- list(call=theCall, threshold=u , map = o,
                    burn = burn, thin = thin, 
                    chains=res ,
                    X.phi = X.phi, X.xi = X.xi,
					acceptance=acc, seed=seed 
                    )
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

  checkEqualsNumeric(cparas, mod.coef, tolerance = tol,msg="gpd: parameter ests page 85 Coles")
  checkEqualsNumeric(cse[2], mod$se[2], tolerance = tol,msg="gpd: standard errors page 85 Coles")
  checkEqualsNumeric(ccov[2,2], mod$cov[2, 2], tolerance = tol,msg="gpd: Covariance page 85 Coles")
  checkEqualsNumeric(cloglik, mod.loglik, tolerance = tol,msg="gpd: loglik page 85 Coles")
  
###################################################################
#   Logical checks on the effect of Gaussian penalization. The smaller the
#    variance, the more the parameter should be drawn towards the
#    mean.

# 2.1 Tests for xi being drawn to 0

  gp1 <- list(c(0, 0), diag(c(10^4, .25)))
  gp2 <- list(c(0, 0), diag(c(10^4, .05)))

  mod1 <- gpd(rain, th=30, priorParameters=gp1)
  mod2 <- gpd(rain, th=30, priorParameters=gp2)

  checkTrue(coef(mod)[2] > coef(mod1)[2],msg="gpd: Gaussian penalization xi being drawn to 0")
  checkTrue(coef(mod1)[2] > coef(mod2)[2],msg="gpd: Gaussian penalization xi being drawn to 0")

# 2.2 Tests for phi being drawn to 0

  gp3 <- list(c(0, 0), diag(c(1, 10^4)))
  gp4 <- list(c(0, 0), diag(c(.1, 10^4)))

  mod3 <- gpd(rain, th=30, priorParameters=gp3)
  mod4 <- gpd(rain, th=30, priorParameters=gp4)

  checkTrue(coef(mod)[1] > coef(mod3)[1],msg="gpd: Gaussian penalization phi being drawn to 0")
  checkTrue(coef(mod3)[1] > coef(mod4)[1],msg="gpd: Gaussian penalization phi being drawn to 0")
  
# 2.3 Tests for xi being drawn to 1
  gp5 <- list(c(0, 1), diag(c(10^4, .25)))
  gp6 <- list(c(0, 1), diag(c(10^4, .05)))

  mod5 <- gpd(rain, th=30, priorParameters=gp5)
  mod6 <- gpd(rain, th=30, priorParameters=gp6)

  checkTrue(1 - coef(mod)[2] > 1 - coef(mod5)[2],msg="gpd: Gaussian penalization xi being drawn to 1")
  checkTrue(1 - coef(mod1)[2] > 1 - coef(mod6)[2],msg="gpd: Gaussian penalization xi being drawn to 1")
  
# 2.4 Tests for phi being drawn to 4 (greater than mle for phi)

  gp7 <- list(c(4, 0), diag(c(1, 10^4)))
  gp8 <- list(c(4, 0), diag(c(.1, 10^4)))

  mod7 <- gpd(rain, th=30, priorParameters=gp7)
  mod8 <- gpd(rain, th=30, priorParameters=gp8)

  checkTrue(4 - coef(mod)[1] > 4 - coef(mod7)[1],msg="gpd: Gaussian penalization phi being drawn to 4")
  checkTrue(4 - coef(mod3)[1] > 4 - coef(mod8)[1],msg="gpd: Gaussian penalization phi being drawn to 4")
  
###################################################################
#   Logical checks on the effect of penalization using lasso or L1 penalization. The smaller the
#    variance, the more the parameter should be drawn towards the
#    mean.

# 2a.1 Tests for xi being drawn to 0

  gp1 <- list(c(0, 0), solve(diag(c(10^4, .25))))
  gp2 <- list(c(0, 0), solve(diag(c(10^4, .05))))

  mod1 <- gpd(rain, th=30, priorParameters=gp1, penalty="lasso")
  mod2 <- gpd(rain, th=30, priorParameters=gp2, penalty="lasso")

  checkTrue(coef(mod)[2] > coef(mod1)[2],msg="gpd: lasso penalization xi being drawn to 0")
  checkTrue(coef(mod1)[2] > coef(mod2)[2],msg="gpd: lasso penalization xi being drawn to 0")

# 2a.2 Tests for phi being drawn to 0

  gp3 <- list(c(0, 0), solve(diag(c(1, 10^4))))
  gp4 <- list(c(0, 0), solve(diag(c(.1, 10^4))))

  mod3 <- gpd(rain, th=30, priorParameters=gp3, penalty="lasso")
  mod4 <- gpd(rain, th=30, priorParameters=gp4, penalty="lasso")

  checkTrue(coef(mod)[1] > coef(mod3)[1],msg="gpd: lasso penalization phi being drawn to 0")
  checkTrue(coef(mod3)[1] > coef(mod4)[1],msg="gpd: lasso penalization phi being drawn to 0")
  
# 2a.3 Tests for xi being drawn to 1
  gp5 <- list(c(0, 1), solve(diag(c(10^4, .25))))
  gp6 <- list(c(0, 1), solve(diag(c(10^4, .05))))

  mod5 <- gpd(rain, th=30, priorParameters=gp5, penalty="lasso")
  mod6 <- gpd(rain, th=30, priorParameters=gp6, penalty="lasso")

  checkTrue(1 - coef(mod)[2] > 1 - coef(mod5)[2],msg="gpd: lasso penalization xi being drawn to 1")
  checkTrue(1 - coef(mod1)[2] > 1 - coef(mod6)[2],msg="gpd: lasso penalization xi being drawn to 1")
  
# 2a.4 Tests for phi being drawn to 4 (greater than mle for phi)

  gp7 <- list(c(4, 0), solve(diag(c(1, 10^4))))
  gp8 <- list(c(4, 0), solve(diag(c(.1, 10^4))))

  mod7 <- gpd(rain, th=30, priorParameters=gp7, penalty="lasso")
  mod8 <- gpd(rain, th=30, priorParameters=gp8, penalty="lasso")

  checkTrue(4 - coef(mod)[1] > 4 - coef(mod7)[1],msg="gpd: lasso penalization phi being drawn to 4")
  checkTrue(4 - coef(mod3)[1] > 4 - coef(mod8)[1],msg="gpd: lasso penalization phi being drawn to 4")
 
########################################################
# Tests on including covariates. Once more, gpd.fit in ismev
# works with sigma inside the optimizer, so we need to tolerate
# some differences and standard errors might be a out.

# 3.0 Reproduce Coles, page 119. Reported log-likelihood is -484.6.

  rtime <- 1:length(rain)
  d <- data.frame(rainfall = rain, time=rtime)
  mod <- gpd(rainfall, th=30, data=d, phi= ~ time, penalty="none")

  checkEqualsNumeric(-484.6, mod$loglik, tolerance = tol,msg="gpd: loglik Coles page 119")
  
####################################################################
# 3.1 Use liver data, compare with ismev. 
#     These are not necessarily sensible models!
#     Start with phi alone.

  mod <- gpd(ALT_M, qu=.7, data=liver,
           phi = ~ ALT_B + dose, xi = ~1,
           penalty="none")

  m <- model.matrix(~ ALT_B + dose, liver)

  ismod <- .ismev.gpd.fit(liver$ALT_M, threshold=quantile(liver$ALT_M, .7), 
                 ydat = m, sigl=2:ncol(m), siglink=exp, show=FALSE)

  checkEqualsNumeric(ismod$mle, coef(mod), tolerance = tol,msg="gpd: covariates in phi only, point ests")
  
# SEs for phi will not be same as for sigma, but we can test xi
  checkEqualsNumeric(ismod$se[length(ismod$se)], mod$se[length(mod$se)], tolerance = tol,msg="gpd: covariates in phi only, standard errors")

######################################################################
# 3.2 Test xi alone.
  mod <- gpd(ALT_M, qu=.7, data=liver,
           phi = ~1, xi = ~ ALT_B + dose,
           penalty="none")

  m <- model.matrix(~ ALT_B + dose, liver)

  ismod <- .ismev.gpd.fit(liver$ALT_M, threshold=quantile(liver$ALT_M, .7), 
                   ydat = m, shl=2:ncol(m), show=FALSE)
  mco <- coef(mod)
  mco[1] <- exp(mco[1])

  checkEqualsNumeric(ismod$mle, mco, tolerance = tol,msg="gpd: covariates in xi only: point ests")
  
# SEs for phi will not be same as for sigma, but we can test xi
  checkEqualsNumeric(ismod$se[-1], mod$se[-1], tolerance = tol,msg="gpd: covariates in xi only: standard errors")

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

  a <- seq(0.1,1,len=10)
  b <- rep(c(-0.5,0.5),each=5)
  data <- makeData(a,b)
  m <- model.matrix(~ a+b, data)
  
  mod <- gpd(y,qu=0.7,data=data,phi=~a,xi=~b,penalty="none")
  ismod <- .ismev.gpd.fit(data$y,threshold=quantile(data$y,0.7),
                   ydat=m,shl=3,sigl=2,siglink=exp, show=FALSE)

  checkEqualsNumeric(ismod$mle, coef(mod),tolerance = tol,msg="gpd: covariates in phi and xi: point ests")
  checkEqualsNumeric(ismod$se, sqrt(diag(mod$cov)),tolerance = tol,msg="gpd: covariates in phi and xi: std errs")

####################################################################
# Check that using priors gives expected behaviour when covariates are included.

# 2.1 Tests for xi being drawn to 0

  b <- rep(c(0.5,1.5),each=5)
  data <- makeData(a=1,b,n=3000)
  
  gp1 <- list(c(0, 0, 0), diag(c(10^4, 0.25, 0.25)))
  gp2 <- list(c(0, 0, 0), diag(c(10^4, 0.25, 0.01)))

  mod0 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,penalty="none")
  mod1 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp1)
  mod2 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp2)

  checkTrue(all(abs(coef(mod0)[2:3]) > abs(coef(mod1)[2:3])),msg="gpd: with covariates, xi drawn to zero")
  checkTrue(abs(coef(mod1)[3]) > abs(coef(mod2)[3]),msg="gpd: with covariates, xi drawn to zero")

# 2.2 Tests for phi being drawn to 0

  a <- seq(0.1,1,len=10)
  data <- makeData(-3 + a,b=-0.1,n=3000)
  data$a <- a
  
  gp4 <- list(c(0, 0, 0), diag(c(1, 1, 10^4)))
  gp5 <- list(c(0, 0, 0), diag(c(0.1, 0.1, 10^4)))

  mod3 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,penalty="none")
  mod4 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp4)
  mod5 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp5)

  checkTrue(all(abs(coef(mod3)[1:2]) > abs(coef(mod4)[1:2])),msg="gpd: with covariates, phi being drawn to 0")
  checkTrue(all(abs(coef(mod4)[1:2]) > abs(coef(mod5)[1:2])),msg="gpd: with covariates, phi being drawn to 0")

# 2.3 Tests for xi being drawn to 2
  b <- rep(c(-0.5,0.5),each=5)
  data <- makeData(a=1,b,n=3000)
  
  gp7 <- list(c(0, 2, 2), diag(c(10^4, 0.25, 0.25)))
  gp8 <- list(c(0, 2, 2), diag(c(10^4, 0.05, 0.05)))

  mod6 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,penalty="none")
  mod7 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp7)
  mod8 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp8)

  checkTrue(all(abs(2 - coef(mod6)[2:3]) > abs(2 - coef(mod7)[2:3])),msg="gpd: with covariates, xi drawn to 2")
  checkTrue(all(abs(2 - coef(mod7)[2:3]) > abs(2 - coef(mod8)[2:3])),msg="gpd: with covariates, xi drawn to 2")

# 2.4 Tests for phi being drawn to 4 

  a <- seq(0.1,1,len=10)
  data <- makeData(2 + a,b=-0.1,n=3000)
  data$a <- a
  
  gp10 <- list(c(0, 4, 0), diag(c(10^4, 1,   10^4)))
  gp11 <- list(c(0, 4, 0), diag(c(10^4, 0.1, 10^4)))

  mod9 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,penalty="none")
  mod10 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp10)
  mod11 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp11)

  checkTrue(abs(4 - coef(mod9)[2])  > abs(4 - coef(mod10)[2]),msg="gpd: with covariates, phi drawn to 4")
  checkTrue(abs(4 - coef(mod10)[2]) > abs(4 - coef(mod11)[2]),msg="gpd: with covariates, phi drawn to 4")

#*************************************************************
  postSum <- function(x){
    t(apply(x$param, 2, function(o){ c(mean=mean(o), se=sd(o)) }))
  } 
  
#************************************************************* 
# 4.1. Test reproducibility

  set.seed(20101110)
  save.seed <- .Random.seed

  set.seed(save.seed)
  bmod <- gpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7), iter=1000,verbose=FALSE, method="sim")
  
  set.seed(save.seed)
  bmod2 <- gpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7), iter=1000,verbose=FALSE, method="sim")
  
  checkEqualsNumeric(bmod$param,bmod2$param,msg="gpd: test simulation reproducibility 1")

  set.seed(bmod$seed)
  bmod3 <- gpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7), iter=1000,verbose=FALSE, method="sim")
  checkEqualsNumeric(bmod$param, bmod3$param,msg="gpd: test simulation reproducibility 2")

#*************************************************************  
# 4.2. Logical test of burn-in

  checkEqualsNumeric(nrow(bmod$chains) - bmod$burn, nrow(bmod$param),msg="gpd: Logical test of burn-in 1")

  iter <- sample(500:1000,1)
  burn <- sample(50,1)
  bmod2 <- gpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7),
                iter=iter, burn=burn, verbose=FALSE, method="sim")

  checkEqualsNumeric(iter-burn, nrow(bmod2$param),msg="gpd: Logical test of burn-in 2")

#*************************************************************
# 4.3. Logical test of thinning

  thin <- 0.5
  iter <- 1000
  bmod <- gpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7),
               iter=iter, thin = thin,verbose=FALSE, method="sim")

  checkEqualsNumeric((nrow(bmod$chains) - bmod$burn) * thin, nrow(bmod$param),msg="gpd: Logical test of thinning 1")

  thin <- 2
  iter <- 1000
  bmod <- gpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7),
               iter=iter, thin = thin, verbose=FALSE, method="sim")

  checkEqualsNumeric((nrow(bmod$chains) - bmod$burn) / thin, nrow(bmod$param),msg="gpd: Logical test of thinning 1")

#*************************************************************  
## Checks of gpd simulation. Centres of distributions and standard deviations 
## should be similar to those coming out of gpd with the same penalty.
## Posterior means are not the same as MAPs and SEs from Hessian are
## approximations, so allow some tolerance.

  tol <- 0.01
  tol.se <- 0.2
# 4.4. Compare MAP and posterior summaries for simple model

  mod <- gpd(ALT_M, data=liver, qu=.7)
  bmod <- gpd(ALT_M, data=liver, th=quantile(liver$ALT_M, .7),verbose=FALSE, method="sim")
  
  checkEqualsNumeric(coef(mod), postSum(bmod)[,1], tolerance=tol,msg="gpd: Compare MAP and posterior summaries for simple model - point ests")
  checkEqualsNumeric(mod$se, postSum(bmod)[,2], tolerance=tol.se,msg="gpd: Compare MAP and posterior summaries for simple model - std errs")

#*************************************************************
# 4.5. Covariates in phi
  tol <- 0.01
  tol.se <- 0.2

  require(MASS,quiet=TRUE) # Use MASS because is ships with R and has no dependencies
              # that don't ship with R.
  liver <- liver
  liver$ndose <- as.numeric(liver$dose)

# function rlm : Fit a linear model by robust regression using an M estimator:

  rmod <- rlm(log(ALT_M) ~ log(ALT_B) + dose, data=liver) 
  liver$resids <- resid(rmod)

  mod <- gpd(resids, data=liver, qu=.7, phi = ~ ndose)

  bmod <- gpd(resids, data=liver, th=quantile(liver$resids, .7),
               phi = ~ ndose,verbose=FALSE, method="sim")
              
  checkEqualsNumeric(coef(mod), postSum(bmod)[,1], tolerance=tol,msg="gpd: Compare MAP and posterior summaries for Covariates in phi - point ests")
  checkEqualsNumeric(mod$se, postSum(bmod)[,2], tolerance=tol.se,msg="gpd: Compare MAP and posterior summaries for Covariates in phi - std errs")

#*************************************************************
# 4.6. Covariates in xi
  tol <- 0.02
  tol.se <- 0.2

  mod <- gpd(resids, data=liver, qu=.7, xi = ~ ndose)
  bmod <- gpd(resids, data=liver, th=quantile(liver$resids, .7),
               xi = ~ ndose,verbose=FALSE, method="sim")

  checkEqualsNumeric(coef(mod), postSum(bmod)[,1], tolerance=tol,msg="gpd: Compare MAP and posterior summaries for Covariates in xi - point ests")
  checkEqualsNumeric(mod$se, postSum(bmod)[,2], tolerance=tol.se,msg="gpd: Compare MAP and posterior summaries for Covariates in xi - std errs")  

#*************************************************************
# 4.7. Covariates in xi and phi
  tol <- 0.02
  tol.se <- 0.3

  mod <- gpd(resids, data=liver, qu=.7,
             xi = ~ ndose, 
             phi = ~ ndose)
  bmod <- gpd(resids, data=liver, th=quantile(liver$resids, .7),
               xi = ~ ndose, 
               phi = ~ ndose,verbose=FALSE, method="sim")
               
  checkEqualsNumeric(coef(mod), postSum(bmod)[,1], tolerance=tol,msg="gpd: Compare MAP and posterior summaries for Covariates in phi and xi - point ests")
  checkEqualsNumeric(mod$se, postSum(bmod)[,2], tolerance=tol.se,msg="gpd: Compare MAP and posterior summaries for Covariates in phi and xi - std errs")  
}

