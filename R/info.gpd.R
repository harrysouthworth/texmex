# TODO: Add comment
# 
# Author: kpzv097
###############################################################################


info.gpd <-
    # Compute the observed information matrix from a gpd object.
    # The expressions are given in Appendix A of Davison & Smith 1990. On page
    # 397 they state s(xTb) = exp(xTb).
	# If penalization is used, the calculation accounts for this, but the resulting
	# estimates of variance will be too low and bias might dominate MSE
function(o, method="observed"){	
    if (class(o) != "gpd"){ stop("object must be of class 'gpd'") }

	if (method != "observed"){
		stop("only 'observed' information is implemented")
	}
	
	# Set out stall
	x <- o$X.phi; z <- o$X.xi
	ns <- ncol(x); nk <- ncol(z)
	s <- coef(o)[1:ns]; k <- -coef(o)[(ns+1):(ns + nk)] # D & S use k = -xi
	s <- exp(colSums(s * t(x))); k <- colSums(k * t(z))
	w <- (o$y - o$threshold) / s

	# Second derivatives of penalties
	p <- matrix(0, nrow=ns+nk, ncol=ns+nk) # OK for MLE, will need a warning for L1
	if (o$penalty %in% c("guassian", "quadratic")){
		Si <- solve(o$priorParameters[[2]])
		for (i in 1:(ns+nk)){
			for (j in 1:(ns + nk)){
				p[i,j] <- 2*Si[i,j]
			}
		}
	} # Close if (o$penatly %in%
	
	# First derivatives of sigma and xi. These are vector derivatives
	d.s <- s * x
	d.k <- z
	
	# Second derivatives of sigma
	d2.s2 <- d.s * x
	# 0 for k because predictor is linear
	
    # First derivatives of loglik. Davison & Smith use k = -xi, so need to adjust
	# derivatives in k by *(-1)
    dl.dk <- (-1 / k^2 * log(1 - k*w) + (1 - 1/k) * w / (1 - k*w) ) * (-1)	
	dl.ds <- ( (w - 1) / (1 - k*w) ) / s
	
	# Second derivatives
	d2l.dk2 <- 2 / k^3 * log(1 - k*w) + 2*w / (k^2 * (1 - k*w)) + (1 - 1/k) * (w / (1 - k*w))^2
	d2l.dkds <- ( (w / s) * (w - 1) / (1 - k*w)^2 ) * (-1)
	d2l.ds2 <- s^(-2) * (1 - 2*w + k*w*w) / (1 - k*w)^2
	
    # Matrix has 4 blocks, 2 of which are the same. Need block for sigma parameters,
	# block for xi parameters and block for the cross of them
	
	Is <- matrix(0, ncol=ns, nrow=ns)
    for (i in 1:ns){
		for (j in 1:ns){ # Lazy coding - only need lower diagonal
			  Is[i,j] <- -sum(x[,i] * x[,j] * (d.s[,i]^2 * d2l.ds2) + d2.s2[,i] * dl.ds)
		}
	}

	Ik <- matrix(0, ncol=nk, nrow=nk)
	for (i in 1:nk){
		for (j in 1:nk){
			Ik[i,j] <- -sum(z[,i] * z[,j] * (d.k[, i]^2 * d2l.dk2 ))
		}
	}
	
	Iks <- matrix(0, ncol=nk, nrow=ns)
	for (i in 1:ns){
		for (j in 1:nk){
			Iks[i,j] <- -sum(z[,j] * x[,i] * d.k[,j] * d.s[,i] * d2l.dkds )
		}
	}

	i <- rbind( cbind(Is, Iks), cbind(t(Iks), Ik))
	res <- try(solve(i - p), silent=TRUE)
	if (class(res) == "try-error"){
	    warning("singular information matrix. Returning a matrix of 0s. Try using bootgpd for inference instead")
	    res <- matrix(0, ncol=ns+nk, nrow=ns+nk)
	}
	res
}

test.info.gpd <- function(){
	lmod <- gpd(log(ALT.M) / log(ALT.B), data=liver, qu=.5, xi=~I(240*as.numeric(dose)), cov="numeric")
	checkTrue(all(sqrt(diag(info.gpd(lmod))) > 0), msg="info.gpd: SDs positive")

	# Check equality to numerical approximation in big samples
	set.seed(20110511)
	tol <- 10^(-3)
	for (i in 1:10){
		x <- rt(10000, 10)
		junk <- gpd(x, qu=.9, penalty="none", cov="numeric")
		msg <- paste("info.gpd: t", i, "equality to numerical", sep="")
		checkEqualsNumeric(junk$cov, info.gpd(junk), tolerance=tol, msg=msg)
	}
}

