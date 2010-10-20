`dgpd` <-
function(x, sigma = 1, xi = 0, u = 0, log=FALSE ){

    n <- length(x)

    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)
    u <- rep(u, length=n)

    if (all(xi == 0)){
        res <- dexp(x-u, sigma, log=TRUE)
    }
    else if (any(xi == 0)){
        res <- numeric(n)
	    wh <- xi == 0
	    res[wh] <- dexp(x[wh] - u[wh], sigma[wh], log=TRUE)
	    res[!wh] <- logb(1 + (xi[!wh] * (x[!wh] - u[!wh]))/sigma[!wh] ) * ( -1/xi[!wh] - 1) - log(sigma[!wh])
    }
    else{
        res <- logb(1 + (xi * (x - u))/sigma ) * ( -1/xi - 1) - log(sigma)
    }

    res <- ifelse( exp(res) >=0, res, 0 )
	
    if (!log){
        res <- exp(res)
    }
    res
}

