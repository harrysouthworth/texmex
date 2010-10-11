`qgpd` <-
function(p , sigma = 1, xi = 1 , u = 0, lower.tail=TRUE, log.p=FALSE ){
    n <- max(length(p), length(sigma), length(xi), length(u))
    p <- rep(p, length=n)
    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)
    u <- rep(u, length=n)

    if (all(xi == 0)){
        res <- qexp(p, 1/sigma)
    }
    else if (any(xi == 0)){
        res <- numeric(length=n)
        wh <- xi == 0
        res[wh] <- qexp(p[wh], 1/sigma[wh])
        res[!wh] <- u[!wh] + (sigma[!wh] / xi[!wh] * p[!wh] ^ (-xi[!wh]) - 1)

    }
	else {
		res <- u + ( sigma * (p^( - xi ) - 1)) / xi
	}
    res
}

