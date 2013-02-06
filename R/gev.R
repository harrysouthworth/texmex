gev <- texmexFamily(name = 'GEV',
                    param = c('mu', 'phi', 'xi'),
                    log.lik = function(data, ...) {
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
                    }, # Close log.lik
                    info = NULL, # will mean that numerical approx gets used
                    delta = function(param, m, model){ # model not used but required by a calling function
                              y <- -log(1 - 1/m)
                              out <- rep(1, 3)

                              out[2] <- -exp(param[2])/param[3] * (1 - y^(-param[3]))
                              out[3] <- exp(param[2]) * param[3]^(-2) * (1 - y^(-param[3])) -
                              exp(param[2])/param[3] * y^(-param[3]) * log(y)
                              out
                    }, # Close delta
                    start = function(data){
                              y <- data$y
                              X.mu <- data$D[[1]]
                              X.phi <- data$D[[2]]
                              X.xi <- data$D[[3]]

                              c(mean(y), rep(.1, ncol(X.mu)-1), log(IQR(y)/2),
                              rep(.01, -1 + ncol(X.phi) + ncol(X.xi)))
                    }, # Close start
                    endpoint = function(param, model){
                                 param[, 1] - exp(param[, 2]) / param[, 3]
                    },
                    rng = function(n, param, model){
                            rgev(n, param[, 1], exp(param[, 2]), param[, 3])
                    },
                    density = function(n, param, model){
                                dgev(n, param[, 1], exp(param[, 2]), param[, 3])
                    },
                    prob = function(n, param, model){
                             pgev(n, param[, 1], exp(param[, 2]), param[, 3])
                    },
                    quant = function(n, param, model){
                              qgev(n, param[, 1], exp(param[, 2]), param[, 3])
                    },
                    resid = function(o){
                              p <- texmexMakeParams(coef(o), o$data$D)
                              scaledY <- (o$data$y - p[, 1]) * p[, 3] / exp(p[, 2])
                              (1 - scaledY)^(-1/xi)
                    }, # Close resid

                    rl = function(m, param, model){ # model not used but required by a calling function
                           param[, 1] - exp(param[, 2])/param[, 3] * (1 - (-log(1 - 1/m))^(-param[, 3]))
                    }
)

