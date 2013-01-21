rgev <- function(n, mu, sigma, xi){
    mu <- rep(mu, length=n)
    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)

    ((-log(runif(n)))^(-xi) - 1)*sigma/xi + mu
}


#x <- rgev(10, 1, 1, .2)
#qgev(pgev(x, 1, 1, .2), 1, 1, .2)
#x



