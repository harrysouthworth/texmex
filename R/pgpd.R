pgpd <-
function(q, sigma = 1, xi = 0, u = 0, lower.tail=TRUE, log.p=FALSE ){

    q <- (q - u) / sigma

    if (xi == 0){
        res <- pexp(q, lower.tail=lower.tail, log.p=log.p )
    }

    else { # Open main else
	    res <- log(1 + xi * q) * (-1/xi) # log survivor function

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
    } # Close main else
	res
}

