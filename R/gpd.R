gpd <-
function(y, data, th, qu, phi= ~1, xi= ~1,
         penalty="gaussian", priorParameters=NULL,
         maxit=10000, trace=0) {

	theCall <- match.call()

    if (class(try(y, silent=TRUE)) == "try-error"){
        if (!missing(data)){
            y <- deparse(substitute(y))
        }
    }

    if (!missing(data)){
        y <- data[, y]
        if (missing(th)){
            th <- quantile(y, qu)
        }
        X.phi <- model.matrix(phi, data)
        X.xi <- model.matrix(xi, data)
        X.phi <- X.phi[y > th,]
        X.xi <- X.xi[y > th,]
    }
    else {
        if (missing(th)){
            th <- quantile(y, qu)
        }
        X.phi <- matrix(ncol=1, rep(1, length(y)))
        X.xi <- matrix(ncol=1, rep(1, length(y)))
        X.phi <- X.phi[y > th, ]
        X.xi <- X.xi[y > th, ]
    }

    rate <- mean(y > th)

    y <- y[y > th]

    if (!is.matrix(X.phi)){
        X.phi <- matrix(X.phi, ncol=1)
    }
    if (!is.matrix(X.xi)){
        X.xi <- matrix(X.xi, ncol=1)
    }

	  if ( length( y ) == 0 ){ stop( "No observations over threshold" ) }
      if ( casefold( penalty ) %in% c( "quadratic" , "gaussian" )  ){
		if ( is.null( priorParameters ) ){
		    priorParameters <- list(rep(0, ncol(X.phi) + ncol(X.xi)),
		                              diag(rep(10^4, ncol(X.phi) + ncol(X.xi)))
		                              )
        }
      }
      else if (casefold(penalty) %in% c("lasso", "l1", "laplace")){
          if (is.null(priorParameters)){
              priorParameters <- list(rep(0, ncol(X.phi) + ncol(X.xi)),
                                      diag(rep(10^(-4), ncol(X.phi) + ncol(X.xi)))
                                      )
          }
          if (length(priorParameters) != 2 | !is.list(priorParameters) ){
              stop("For Laplace prior, priorParameters should be a list
                    of length 2, the second element of which should be
                    a diagonal matrix")
          }
          if (!is.matrix(priorParameters[[2]])){
              # Lasso penalty is a bit unintuitive. Try to help out.
              priorParameters[[2]] <- diag(rep(priorParameters[[2]], ncol(X.phi) + ncol(X.xi)))      
          }
          if (!all(wh == diag(diag(wh)))){
              warning("some off-diagonal elements of the covariance are non-zero. Only
                       the diagonal is used in penalization")
          }
      }

        if (!is.null(priorParameters)){
            dimp <- ncol(X.phi) + ncol(X.xi)
            if (length(priorParameters[[1]]) != dimp){
                stop("wrong number of parameters in prior (doesn't match phi and xi formulas)")
            }
            else if (length(diag(priorParameters[[2]])) != dimp){
                stop("wrong dimension of prior covariance (doesn't match phi and xi formulas)")
            }

        }

        gpd.lik <- function(par, y, th, X.phi, X.xi,
                            penalty = "none", priorParameters=NULL,
                            maxit=maxit, trace=trace) {
            keepsc <- par[1:ncol(X.phi)]
            keepxi <- par[ - (1:ncol(X.phi))]

            sc <- colSums(par[1:ncol(X.phi)] * t(X.phi))
            xi <- colSums(par[(ncol(X.phi) + 1):(ncol(X.phi) + ncol(X.xi))] * t(X.xi))

            y <- (y - th)/exp(sc)
            y <- 1 + xi * y

            if (min(y) <= 0) 
                l <- 10^6
            else l <- sum(sc) + sum(log(y) * (1/xi + 1)) # neg log lik
            if (casefold(penalty) == "none") 
                l <- l
            else if (casefold(penalty) %in% c("quadratic", "gaussian")) {

                p <- mahalanobis(matrix(c(keepsc, keepxi), nrow = 1), 
                  center = priorParameters[[1]], cov = priorParameters[[2]])
                l <- l + p
            }
            else if (casefold(penalty) %in% c("lasso", "l1", "laplace")){
                p <- sum(abs(c(keepsc, keepxi) - priorParameters[[1]]) * diag(priorParameters[[2]]))
                l <- l + p
            }
            else stop("penalty can be 'none', 'lasso' or 'gaussian'")
            l
        }

        start <- c(log(mean(y)), rep(.00001, -1 + ncol(X.phi) + ncol(X.xi)))

        o <- optim(par = start, fn = gpd.lik,
                   y = y, X.phi = X.phi, X.xi = X.xi, th = th, 
                   penalty = penalty, 
                   control = list(maxit = maxit, trace = trace), 
                   priorParameters = priorParameters,
		hessian =  TRUE )

      if (o$convergence != 0){
          warning("Non-convergence in gpd")
      }
      o$cov <- solve(o$hessian)
	  o$hessian <- NULL
	  o$se <- sqrt( diag( o$cov ) )
	  o$threshold <- th
	  o$penalty <- penalty
	  o$mle <- o$par
      names(o$mle) <- c(paste("phi:", colnames(X.phi)),
                        paste("xi:", colnames(X.xi)))
	  o$par <- NULL
	  o$rate <- rate
	  o$call <- theCall
	  o$y <- y
	  o$X.phi <- X.phi
	  o$X.xi <- X.xi
	  if (missing(data)) { data <- NULL }
      o$data <- data
	  o$loglik <- -gpd.lik(o$mle, y, th, X.phi, X.xi)
    o$value <- NULL # since this duplicates the $loglik value
    o$counts <- NULL # not required
	  oldClass( o ) <- "gpd"
	  o
    }

test(gpd) <- function(){
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
  
  mod.coef[1] <- exp(mod.coef)
  names(mod.coef)[1] <- "sigma"
  
  mod.loglik <- mod$loglik
  mod.cov22 <- mod$cov[2, 2]

  checkEqualsNumeric(mod.coef, cparas, tolerance = tol)
  checkEqualsNumeric(mod$se[2], cse[2], tolerance = tol)
  checkEqualsNumeric(mod$cov[2, 2], ccov[2,2], tolerance = tol)
  checkEqualsNumeric(mod.loglik, cloglik, tolerance = tol)
  
###################################################################
#   Logical checks on the effect of penalization. The smaller the
#    variance, the more the parameter should be drawn towards the
#    mean.

# 2.1 Tests for xi being drawn to 0

  gp1 <- list(c(0, 0), diag(c(10^4, .25)))
  gp2 <- list(c(0, 0), diag(c(10^4, .05)))

  mod1 <- gpd(rain, th=30, priorParameters=gp1)
  mod2 <- gpd(rain, th=30, priorParameters=gp2)

  checkTrue(coef(mod)[2] > coef(mod1)[2])
  checkTrue(coef(mod1)[2] > coef(mod2)[2])

# 2.2 Tests for phi being drawn to 0

  gp3 <- list(c(0, 0), diag(c(1, 10^4)))
  gp4 <- list(c(0, 0), diag(c(.1, 10^4)))

  mod3 <- gpd(rain, th=30, priorParameters=gp3)
  mod4 <- gpd(rain, th=30, priorParameters=gp4)

  checkTrue(coef(mod)[1] > coef(mod3)[1])
  checkTrue(coef(mod3)[1] > coef(mod4)[1])
  
# 2.3 Tests for xi being drawn to 1
  gp5 <- list(c(0, 1), diag(c(10^4, .25)))
  gp6 <- list(c(0, 1), diag(c(10^4, .05)))

  mod5 <- gpd(rain, th=30, priorParameters=gp5)
  mod6 <- gpd(rain, th=30, priorParameters=gp6)

  checkTrue(1 - coef(mod)[2] > 1 - coef(mod5)[2])
  checkTrue(1 - coef(mod1)[2] > 1 - coef(mod6)[2])
  
# 2.4 Tests for phi being drawn to 4 (greater than mle for phi)

  gp7 <- list(c(4, 0), diag(c(1, 10^4)))
  gp8 <- list(c(4, 0), diag(c(.1, 10^4)))

  mod7 <- gpd(rain, th=30, priorParameters=gp7)
  mod8 <- gpd(rain, th=30, priorParameters=gp8)

  checkTrue(4 - coef(mod)[1] > 4 - coef(mod7)[1])
  checkTrue(4 - coef(mod3)[1] > 4 - coef(mod8)[1])  
  
########################################################
# Tests on including covariates. Once more, gpd.fit in ismev
# works with sigma inside the optimizer, so we need to tolerate
# some differences and standard errors might be a out.

# 3.0 Reproduce Coles, page 119. Reported log-likelihood is -484.6.

  rtime <- 1:length(rain)
  d <- data.frame(rainfall = rain, time=rtime)
  mod <- gpd(rainfall, th=30, data=d, phi= ~ time, penalty="none")

  checkEqualsNumeric(mod$loglik, -484.6, tolerance = tol)
  
####################################################################
# 3.1 Use liver data, compare with ismev. 
#     These are not necessarily sensible models!
#     Start with phi alone.

  mod <- gpd(ALT_M, qu=.7, data=liver,
           phi = ~ ALT_B + dose, xi = ~1,
           penalty="none")

  m <- model.matrix(~ ALT_B + dose, liver)

  ismod <- gpd.fit(liver$ALT_M, th=quantile(liver$ALT_M, .7), 
                 ydat = m, sigl=2:ncol(m), siglink=exp, show=FALSE)

  checkEqualsNumeric(coef(mod), ismod$mle, tolerance = tol)
  
# SEs for phi will not be same as for sigma, but we can test xi
  checkEqualsNumeric(ismod$se[length(ismod$se)], mod$se[length(mod$se)], tolerance = tol)

######################################################################
# 3.2 Test xi alone.
  mod <- gpd(ALT_M, qu=.7, data=liver,
           phi = ~1, xi = ~ ALT_B + dose,
           penalty="none")

  m <- model.matrix(~ ALT_B + dose, liver)

  ismod <- gpd.fit(liver$ALT_M, th=quantile(liver$ALT_M, .7), 
                   ydat = m, shl=2:ncol(m), show=FALSE)
  mco <- coef(mod)
  mco[1] <- exp(mco[1])

  checkEqualsNumeric(mco, ismod$mle, tolerance = tol)
  
# SEs for phi will not be same as for sigma, but we can test xi
  checkEqualsNumeric(ismod$se[-1], mod$se[-1], tolerance = tol)

######################################################################
# 3.3 Test phi & xi simultaneously. Use simulated data.

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
  ismod <- gpd.fit(data$y,th=quantile(data$y,0.7),
                   ydat=m,shl=3,sigl=2,siglink=exp, show=FALSE)

  checkEqualsNumeric(coef(mod), ismod$mle, "Test phi & xi simultaneously, coefs\n",tolerance = tol)
  checkEqualsNumeric(sqrt(diag(mod$cov)), ismod$se, "Test phi & xi simultaneously, se\n",tolerance = tol)

####################################################################
# Check that using priors gives expected behaviour when covariates are included.

# 2.1 Tests for xi being drawn to 0

  b <- rep(c(0.5,1.5),each=5)
  data <- makeData(a=1,b,n=3000)
  
  gp1 <- list(c(0, 0, 0), diag(c(10^4, 0.25, 0.25)))
  gp2 <- list(c(0, 0, 0), diag(c(10^4, 0.05, 0.05)))

  mod0 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,penalty="none")
  mod1 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp1)
  mod2 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp2)

  checkTrue(abs(coef(mod0)[2:3]) > abs(coef(mod1)[2:3]))
  checkTrue(abs(coef(mod1)[2:3]) > abs(coef(mod2)[2:3]))

# 2.2 Tests for phi being drawn to 0

  a <- seq(0.1,1,len=10)
  data <- makeData(-3 + a,b=-0.1,n=3000)
  data$a <- a
  
  gp4 <- list(c(0, 0, 0), diag(c(1, 1, 10^4)))
  gp5 <- list(c(0, 0, 0), diag(c(0.1, 0.1, 10^4)))

  mod3 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,penalty="none")
  mod4 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp4)
  mod5 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp5)

  checkTrue(abs(coef(mod3)[1:2]) > abs(coef(mod4)[1:2]))
  checkTrue(abs(coef(mod4)[1:2]) > abs(coef(mod5)[1:2]))

# 2.3 Tests for xi being drawn to 2
  b <- rep(c(-0.5,0.5),each=5)
  data <- makeData(a=1,b,n=3000)
  
  gp7 <- list(c(0, 2, 2), diag(c(10^4, 0.25, 0.25)))
  gp8 <- list(c(0, 2, 2), diag(c(10^4, 0.05, 0.05)))

  mod6 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,penalty="none")
  mod7 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp7)
  mod8 <- gpd(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp8)

  checkTrue(abs(2 - coef(mod6)[2:3]) > abs(2 - coef(mod7)[2:3]))
  checkTrue(abs(2 - coef(mod7)[2:3]) > abs(2 - coef(mod8)[2:3]))

# 2.4 Tests for phi being drawn to 4 

  a <- seq(0.1,1,len=10)
  data <- makeData(2 + a,b=-0.1,n=3000)
  data$a <- a
  
  gp10 <- list(c(0, 4, 0), diag(c(10^4, 1,   10^4)))
  gp11 <- list(c(0, 4, 0), diag(c(10^4, 0.1, 10^4)))

  mod9 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,penalty="none")
  mod10 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp10)
  mod11 <- gpd(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp11)

  checkTrue(abs(4 - coef(mod9)[2])  > abs(4 - coef(mod10)[2]))
  checkTrue(abs(4 - coef(mod10)[2]) > abs(4 - coef(mod11)[2]))

}