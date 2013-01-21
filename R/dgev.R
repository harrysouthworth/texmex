dgev <- function(x, sigma, xi, log.d=FALSE){
    n <- length(x)
    
    mu <- rep(mu, length=n)
    sigma <- rep(sigma, length=n)
    xi <- rep(xi, length=n)

    a <- function(x, m, s, z){
        1 + z / s * (x - m)
    }

    a <- a(x, mu, sigma, xi)

    ll <- -n*log(sigma) - (1 + 1/xi) * (sum(log(a)) - sum(a^(-1/xi)))
    if (log.d){ ll }
    else { exp(ll) }
}
