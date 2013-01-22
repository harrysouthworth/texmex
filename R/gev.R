gev.loglik <- function(data, ...) {
  y <- data$y
  X.mu <- data$D$mu
  X.phi <- data$D$phi
  X.xi <- data$D$xi

  n.mu <- ncol(X.mu)
  n.phi <- n.mu + ncol(X.phi)
  n.end <- n.phi + ncol(X.xi)

  function(param) {
    stopifnot(length(param) == n.end)
    mu <- X.mu %*% param[1:n.mu]
    phi <- X.phi %*% param[(1 + n.mu):n.phi]
    xi <- X.xi %*% param[(1 + n.phi):n.end]
    sum(dgev(y, mu, exp(phi), xi, log.d=TRUE))
  }
}

rl.gev <- function(m, param, model){ # model not used but required by a calling function
    param[1] - exp(param[2])/param[3] * (1 - (-log(1 - 1/m))^(-param[3]))
}

gev.delta <- function(param, m, model){ # model not used but required by a calling function
    y <- -log(1 - 1/m)
    out <- rep(1, 3)

    out[2] <- -1/param[3] * (1 - y^(-xi))
    out[3] <- param[2] * param[3]^(-2) * (1 - y^(-xi)) - param[2]/param[3] * y^(-xi) * log(y)
    out
}

gev.start <- function(data){
    y <- data$y
    X.mu <- data$D[[1]]
    X.phi <- data$D[[2]]
    X.xi <- data$D[[3]]

    c(mean(y), rep(.1, ncol(X.mu)-1), log(IQR(y)/2),
      rep(.01, -1 + ncol(X.phi) + ncol(X.xi)))
}

gev.rng <- function(n, param, model){
   rgev(n, param[1], exp(param[2]), param[3])
}
gev.dens <- function(n, param, model){
   dgev(n, param[1], exp(param[2]), param[3])
}
gev.prob <- function(n, param, model){
    pgev(n, param[1], exp(param[2]), param[3])
}
gev.quant <- function(n, param, model){
    qgev(n, param[1], exp(param[2]), param[3])
}

gev <- list(name = 'GEV',
            log.lik = gev.loglik,
            param = c('mu', 'phi', 'xi'),
            info = NULL,
            start = gev.start,
            rng = gev.rng,
            density = gev.dens,
            prob = gev.prob,
            quant = gev.quant,
            resid=function(o){ NULL },
            rl=rl.gev)

