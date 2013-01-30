texmexFamily <-
    # Create an object of class 'texmexFamily'. It is allowable to have
    # info, start and resid as NULL, but all other information must
    # be provided by the user.
function(name, log.lik, param, info=NULL, start=NULL, resid=NULL,
                         rl, delta, density, rng, prob, quant){
    res <- list(name=name, log.lik=log.lik, param=param, info=info, start=start,
                resid=resid, rl=rl, delta=delta, density=density, rng=rng,
                prob=prob, quant=quant)

    oldClass(res) <- 'texmexFamily'
    res
}


print.texmexFamily <- function(x, verbose=TRUE, ...){
    if (is.null(x$info)){ info <- 'Numerical approximation' }
    else { info <- 'Closed form' }

    if (is.null(x$start)){ start <- 'Random' }
    else { start <- 'Data dependent' }

    cat('Family:      ', x$name, '\n')
    if (verbose){
        cat('Parameters:  ', x$param, '\n')
        cat('Information: ', info, '\n')
    }

    invisible()
}
