info.gpd <-
    # Compute the observed information matrix from a gpd object.
    # The expressions are given in Appendix A of Davison & Smith 1990.
    # Note we are using a simpler parameterisation in which phi = log(sigma)
    # and xi are both linear in their covariates. Xi is -k used in Davison and Smith.
	# If penalization is used, the calculation accounts for this, but the resulting
	# estimates of variance will be too low and bias might dominate MSE
function(o, method="observed"){
  if (class(o) != "gpd"){ stop("object must be of class 'gpd'") }
	if (method != "observed"){ stop("only 'observed' information is implemented") }

	x <- o$X.phi; z <- o$X.xi
	ns <- ncol(x); nk <- ncol(z)
	phi <- coef(o)[1:ns]
  xi <- coef(o)[(ns+1):(ns + nk)]

	phi.i <- colSums(phi * t(x))
  xi.i <- colSums(xi * t(z))
	w.i <- (o$y - o$threshold) / exp(phi.i)

# Second derivatives of penalties
	p <- matrix(0, nrow=ns+nk, ncol=ns+nk)
	if (o$penalty %in% c("gaussian", "quadratic")){ # note if Lasso penalty used then 2nd deriv is zero hence no term for this
		Si <- solve(o$priorParameters[[2]])
		for (i in 1:(ns+nk)){
			for (j in 1:(ns + nk)){
				p[i,j] <- 2*Si[i,j]
			}
		}
	}

# Second and mixed derivatives of log-lik wrt coefficients of linear predictors

  d2li.dphi2 <- -(1 + 1/xi.i) * xi.i * w.i / (1 + xi.i*w.i)^2
  d2li.dphidxi <- 1/xi.i^2 * (1/(1 + xi.i*w.i) - 1) + (1+1/xi.i)*w.i/(1 + xi.i*w.i)^2
  d2li.dxi2 <- -2/xi.i^3 * log(1 + xi.i*w.i) + 2*w.i/(xi.i^2 * (1 + xi.i*w.i)) + (1 + 1/xi.i)*w.i^2/(1 + xi.i*w.i)^2

# Matrix has 4 blocks, 2 of which are transposes of each other. Need block for phi parameters,
# block for xi parameters and block for the cross of them.

  Ip <- matrix(0, ncol=ns, nrow=ns)
  for (u in 1:ns){
	  for (v in 1:ns){
		  Ip[u,v] <- -sum(x[,u] * x[,v] * d2li.dphi2)
	  }
	}

	Ix <- matrix(0, ncol=nk, nrow=nk)
	for (s in 1:nk){
		for (t in 1:nk){
			Ix[s,t] <- -sum(z[,s] * z[,t] * d2li.dxi2)
		}
	}

	Ipx <- matrix(0, ncol=nk, nrow=ns)
	for (u in 1:ns){
		for (s in 1:nk){
			Ipx[u,s] <- -sum(z[,s] * x[,u] * d2li.dphidxi )
		}
	}

	i <- rbind( cbind(Ip, Ipx), cbind(t(Ipx), Ix))

	# return observed Information matrix.   Note that an estimate of the covariance matrix is given by the inverse of this matrix.
  i - p
}

test.info.gpd <- function(){
	lmod <- gpd(ALT.M, data=liver, qu=.5, xi=~I(240*as.numeric(dose)), cov="numeric")
	checkTrue(all(sqrt(diag(solve(info.gpd(lmod)))) > 0), msg="info.gpd: SDs positive")

	# Check equality to numerical approximation in big samples
	set.seed(20110923)
	tol <- 10^(-3)
	for (i in 1:10){
            x <- rt(10000, 10)
            junk <- gpd(x, qu=.9, penalty="none", cov="numeric")
            msg <- paste("info.gpd: t", i, "equality to numerical", sep="")
            checkEqualsNumeric(junk$cov, solve(info.gpd(junk)), tolerance=tol, msg=msg)

  # check estimation when we have a penalty
    gp1 <- list(c(0, 0), diag(c(10^4, .05)))
    gp2 <- list(c(0, 0), diag(c(.1, 10^4)))
		junk1 <- gpd(x, qu=.9, priorParameters = gp1, cov="numeric")
		junk2 <- gpd(x, qu=.9, priorParameters = gp2, cov="numeric")
		msg1 <- paste("info.gpd: t", i, "equality to numerical, penalty on xi", sep="")
		msg2 <- paste("info.gpd: t", i, "equality to numerical, penalty on phi", sep="")
    tol <- 0.01
		checkEqualsNumeric(junk1$cov, solve(info.gpd(junk1)), tolerance=tol, msg=msg1)
		checkEqualsNumeric(junk2$cov, solve(info.gpd(junk2)), tolerance=tol, msg=msg2)

  # check estimation when we have covariates
    n <- 10000
    x <- 1/runif(n)
    data <- data.frame(x=x,y=rexp(n,exp(2 + x)))

    junk3 <- gpd(y,data=data,phi =~ x,th=0)
    msg3 <- paste("info.gpd: t",i,"equality to numerical, covariates in phi",sep="")
    checkEqualsNumeric(junk3$cov, solve(info.gpd(junk3)), tolerance=tol, msg=msg3)

    x <- runif(n,-0.5,0.5)
    data <- data.frame(x=x,y = rgpd(n,sigma = exp(3+2*x), xi=x))

    junk4 <- gpd(y,data=data,phi=~x, xi = ~ x,th=0)
    msg4 <- paste("info.gpd: t",i,"equality to numerical, covariates in phi and xi",sep="")
    checkEqualsNumeric(junk4$cov, solve(info.gpd(junk4)), tolerance=tol, msg=msg4)
  }
}

