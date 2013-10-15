info.gev <-
  # Compute the observed information matrix for a gev family.
  # The form of the information matrix is given by Prescott and Walden, 1980.
  # Note we are using a parameterisation in which phi = log(alpha)
  # and xi are both linear in their covariates. Xi is -k used in Prescott and Walden.
  # If penalization is used, the calculation accounts for this, but the resulting
  # estimates of variance will be too low and bias might dominate MSE
function(o, method="observed"){
    if (class(o) != "gpd"){ stop("object must be of class 'gpd'") }
    if (method != "observed"){ stop("only 'observed' information is implemented") }

    w <- o$data$D$mu; x <- o$data$D$phi; z <- o$data$D$xi
    nm <- ncol(w); ns <- ncol(x); nk <- ncol(z)

    mu <- coef(o)[1:nm]
    phi <- coef(o)[(nm + 1):(nm + ns)]
    xi <- coef(o)[(nm + ns + 1):(nm + ns + nk)]

    mu.i <- colSums(mu * t(w))
    phi.i <- colSums(phi * t(x))
    xi.i <- colSums(xi * t(z))
    w.i <- (o$data$y - o$threshold) / exp(phi.i)

    # Warn if regularity condition fails
    if (any(xi.i < -.50)){ warning("Fitted values of xi < -0.5") }

    p <- (1 + xi.i)^2 * gamma(1 + 2*xi.i)
    q <- gamma(2 + xi.i) * (digamma(1 + xi.i) + (1 + xi.i)/xi.i)
    g <- -digamma(1) # \gamma in Prescott & Walden

    # Get second derivatives of penalty terms
    pen <- matrix(0, nrow=nm+ns+nk, ncol=nm+ns+nk)
    if (o$penalty %in% c("gaussian", "quadratic")){ # note if Lasso penalty used then 2nd deriv is zero hence no term for this
      Si <- solve(o$priorParameters[[2]])
      for (i in 1:(nm+ns+nk)){
        for (j in 1:(nm+ns+nk)){
          p[i,j] <- 2*Si[i,j]
        }
      } # Close for i
    } # Close if

    # Second derivatives
    d2li.dmu2 <- p/(exp(phi.i)^2)
    d2li.dxi2 <- 1/(xi.i*xi.i) * (pi*pi / 6 + (1 - g + 1/xi.i)^2 - 2*q/xi.i + p/(xi.i*xi.i))
    d2li.dphi2 <- 1/(exp(2*phi.i) * xi.i^2) * (1 - 2*gamma(2 + xi.i) + p)
    
    d2li.dmudphi <- -1/(exp(2*phi.i) * xi.i) * (p - gamma(2 + xi.i))
    d2li.dmudxi <- 1/(exp(phi.i) * xi.i) * (q - p/xi.i)
    d2li.dphidxi <- 1 / (exp(phi.i) * xi.i^2) * (1 - g + (1 - gamma(2 + xi.i))/xi.i - q + p/xi.i)





}