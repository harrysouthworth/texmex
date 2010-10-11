`rgpd` <-
function(n, sigma = 1, xi = 0 , u = 0 ){

	# Check parameter vectors have correct length
    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)
    u <- rep(u, length=n)

    # Get random numbers
    if ( all( xi == 0 ) ){
		res <- rexp( n, 1/sigma )+ u
	}
	else if (any(xi == 0)){
	    res <- numeric(n)

	    wh <- xi == 0
	    res[wh] <- rexp(sum(xi == 0), 1/sigma[wh]) + u[wh]

	    res[!wh] <- sigma[!wh]/xi[!wh] * (runif(sum(!wh))^(-xi[!wh]) - 1) + u[!wh]
	}
	else {
		res <- u + sigma / xi  * (runif(n)^( -xi ) - 1)
    }
    res
}

