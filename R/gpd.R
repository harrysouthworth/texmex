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

gpd <- function(){
    res <-  list(name = 'GPD',
                 log.lik = gpd.loglik,
                 param = c('phi', 'xi'),
                 info = gpd.info,
                 start = gpd.start,
                 resid = gpd.residuals,
                 density=dgpd,
                 rng=rgpd)
    oldClass(res) <- 'texmexFamily'
    res
}
