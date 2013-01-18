gpd.loglik <- function(data, th, ...) {
  y <- data$y
  X.phi <- data$D$phi
  X.xi <- data$D$xi

  n.phi <- ncol(X.phi)
  n.end <- n.phi + ncol(X.xi)
  function(param) {
    stopifnot(length(param) == n.end)
    phi <- X.phi %*% param[1:n.phi]
    xi <- X.xi %*% param[(1 + n.phi):n.end]
    sum(dgpd(y, exp(phi), xi, u=th, log.d=TRUE))
  }
}

gpd.start <- function(data){
    y <- data$y
    X.phi <- data$D[[1]]
    X.xi <- data$D[[2]]

    c(log(mean(y)), rep(1e-05, -1 + ncol(X.phi) + ncol(X.xi)))
}

gpd.residuals <- function(o){
    fittedScale <- fittedGPDscale(o)
    fittedShape <- fittedGPDshape(o)
    scaledY <- fittedShape * (o$data$y - o$threshold) / fittedScale
    c(1/fittedShape * log(1 + scaledY)) # Standard exponential
}

gpd.delta <- function(A, K){
   # This is not exact if a prior (penalty) function is used, but
   # the CI is approximate anyway.

    out <- matrix(0, nrow=2, ncol=length(K))

    if (A[3] == 0){ # exponential case
        out[1,] <- exp(A[2]) * log(K * A[1])
    } else {
        out[1,] <- exp(A[2]) / A[3] * ((K*A[1])^A[3] - 1)
        out[2,] <- -exp(A[2]) / (A[3]*A[3]) * ( (K * A[1] )^A[3] - 1 ) +
                   exp(A[2]) / A[3] * (K * A[1])^A[3] * log(K * A[1])
    }

   out
}


gpd <- list(name = 'GPD',
            log.lik = gpd.loglik,
            param = c('phi', 'xi'),
            info = gpd.info,
            start = gpd.start,
            resid = gpd.residuals,
            delta = gpd.delta,
            density=dgpd,
            rng=rgpd,
            prob=pgpd,
            quant=qgpd)
oldClass(gpd) <- 'texmexFamily'


