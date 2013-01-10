gpd.loglik <- function(data) {
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




gpd <- function(){
  res <-  list(name = 'GPD',
               log.lik = gpd.loglik,
               param = c('phi', 'xi'),
               info = gpd.info)
  oldClass(res) <- 'texmexFamily'
  res
}
                                  #,             # function(o, ...)
#            start = gpd.start            # function
#            resid = gpd.resid,           # function       )
