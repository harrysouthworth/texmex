rgev <- function(n, mu, sigma, xi){

    mu <- rep(mu, length=n)
    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)

    if( xi[1] != 0 ) {
      res <- ((-log(runif(n)))^(-xi) - 1)*sigma/xi + mu
    } else {
      res <- mu - sigma*(log(-log(runif(n))))
    }
    res
}


#x <- rgev(10, 1, 1, .2)
#qgev(pgev(x, 1, 1, .2), 1, 1, .2)
#x



