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

gev.start <- NULL

