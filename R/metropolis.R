texmexMetropolis <- function(x, log.lik, proposals, verbose, trace){
    last.cost <- log.lik(x[1,])
    if (!is.finite(last.cost)) {
      stop("Start is infeasible.")
    }

    acc <- 0
    for(i in 2:nrow(x)){
      if( verbose){
        if(i %% trace == 0) cat(i, " steps taken\n" )
      }
      prop <- proposals[i - 1,] + x[i - 1,]
      top <- log.lik(prop)
      delta <- top - last.cost
      if (is.finite(top) && ((delta >= 0) ||
                             (runif(1) <= exp(delta)))) {
        x[i, ] <- prop
        last.cost <- top
        acc <- 1 + acc
      }
      else {
        x[i, ] <- x[i-1,]
      }
    } # Close for(i in 2:nrow

    acc <- acc / nrow(x)
    if (acc < .1) {
        warning("Acceptance rate in Metropolis algorithm is low.")
    }
    if ((trace < nrow(x)) & verbose) {
        cat("Acceptance rate:", round(acc , 3) , "\n")
    }

    attr(x, "acceptance") <- acc
    x
}
