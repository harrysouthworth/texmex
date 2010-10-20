pgpd <-
function(q, sigma = 1, xi = 0, u = 0, lower.tail=TRUE, log.p=FALSE ){

    q <- (q - u) / sigma
    n <- length(q)

    xi <- rep(xi, n)
    sigma <- rep(sigma, n)
    u <- rep(u, n)

    if (all(xi == 0)){
        res <- pexp(q, log.p=TRUE, lower.tail = FALSE)
    }
    else if (any(xi == 0)){
        res <- numeric(n)
        wh <- xi == 0
        res[wh] <- pexp(q[wh], log.p=TRUE , lower.tail=FALSE)
        res[!wh] <- log(1 + xi[!wh]*q[!wh]) * (-1/xi[!wh])
    }
    else {
        res <- log(1 + xi * q) * (-1/xi)
    }

    if (!log.p){
        res <- exp(res) # survivor function
        if (lower.tail){
            res <- 1 - res
        }
    } # Close if (!log.p
    else { # if want log
        if(lower.tail){
            res <- log(1 - exp(res))
        }
    }

	res
}

