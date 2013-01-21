rgev <- function(n, mu, sigma, xi){
    mu <- rep(mu, length=n)
    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)

    ((-log(runif(n)))^(-xi) - 1)*sigma/xi - mu
}
