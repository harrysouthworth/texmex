gpd.loglik <- function(data, ...) {
  dots <- as.list(substitute(list(...)))[-1]
  if (!is.element('th', names(dots))) {
    th <- quantile(data$y, dots$qu)
  }

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


gpd <- function(){
  res <-  list(name = 'GPD',
               log.lik = gpd.loglik,
               param = c('phi', 'xi'),
               info = gpd.info,
               start = gpd.start)
  oldClass(res) <- 'texmexFamily'
  res
}

#            resid = gpd.resid,           # function       )

print.texmexFamily <- function(x, ...){
    cat('Family: ', x$name, '\n')
    cat('Param: ', x$param)
    invisible()
}
