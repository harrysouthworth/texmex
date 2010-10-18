`dgpd` <-
function(x, sigma = 1, xi = 0, u = 0, log=FALSE ){
	if ( xi != 0 ){
		res <- logb(1 + (xi * (x - u))/sigma ) * ( -1/xi - 1) - log(sigma)
		res <- ifelse( exp(res) >=0, res, 0 )
	}
	else {
	    res <- dexp( x-u, sigma, log=log )
	}
    if (!log){
        res <- exp(res)
    }
    res
}

