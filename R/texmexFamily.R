texmexFamily <- function(name, log.lik, param,
                         info=NULL, start=NULL, resid){
    res <- function(){
        if (is.null(info)) { info <- function(){ NULL }}
        if (is.null(start)) { start <- function(){ NULL }}

        o <- list(name=name, log.lik=log.lik, param=param, info=info,
                  start=start, resid=resid)
        oldClass(o) <- 'texmexFamily'
        o
    }
    res
}


print.texmexFamily <- function(x, ...){
    if (is.null(x$info)){ info <- 'Numerical approximation' }
    else { info <- 'Closed form' }

    if (is.null(x$start)){ start <- 'Random' }
    else { start <- 'Data dependent' }

    cat('Family:      ', x$name, '\n')
    cat('Parameters:  ', x$param, '\n')
    cat('Information: ', info, '\n')
    cat('Start:       ', start, '\n')

    
    invisible()
}
