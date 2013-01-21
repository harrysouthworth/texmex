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

# delta, rl, resid

rl.gev <- function(m, param, model){ # model not used but required by a calling function
    param[1] - param[2]/param[3] * (1 - (-log(1 - 1/m))^(-1/param[3]))
}

gev.delta <- function(param, m, model){ # model not used but required by a calling function
    y <- -log(1 - 1/m)
    out <- rep(1, 3)

    out[2] <- -1/param[3] * (1 - y^(-xi))
    out[3] <- param[2] * param[3]^(-2) * (1 - y^(-xi)) - param[2]/param[3] * y^(-xi) * log(y)
    out
}

gev <- list(name = 'GEV',
            loglik = gev.loglik,
            param = c('mu', 'phi', 'xi'),
            info = NULL,
            start = NULL,
            rng = rgev,
            density = dgev,
            prob = pgev,
            quant = qgev,
            resid=NULL)
