#' @include texmexFamily.R
#' @export gumbel
NULL

gumbel <- texmexFamily(name = 'Gumbel',
					log.lik = function(data, ...) {
						y <- data$y
						X.mu <- data$D$mu 
						X.phi <- data$D$phi #  phi = log(sigma)
						
						n.mu <- ncol(X.mu)
						n.end <- n.mu + ncol(X.phi)
						
						function(param) {
							stopifnot(length(param) == n.end)
							mu <- X.mu %*% param[1:n.mu]
							phi <- X.phi %*% param[(1 + n.mu):n.end]
							z <- (y-mu)/exp(phi)
							sum(-exp(-z)-z-phi)
						}
					}, # Close log.lik
					param = c(mu=0, phi=0),
					info = NULL,
					sandwich = NULL,
					start = function(data){
						y <- data$y
						X.mu <- data$D$mu
						X.phi <- data$D$phi
						c(mean(y), rep(0, ncol(X.mu) - 1), 
						  log(IQR(y)/2), rep(0.001, ncol(X.phi)-1))
					}, # Close start
					
					resid = function(o){
						p <- texmexMakeParams(coef(o), o$data$D)
						(o$data$y - p[,1]) / exp(p[,2])  # Standard gumbel
					}, # Close resid
					
					endpoint = function(param, model){
						Inf
					},
					rl = function(m, param, model){
						param[,1] - exp(param[,2])* log(-log(1-1/m))
					}, # close rl
					delta = function(param,m,model){ # follows argument in Coles eqn (4.15) for GPD
						out <- rep(1,2) # can have vector output as call only ever assumes m length 1

						out[1] <- 1
						out[2] <- -exp(param[2]) * log(-log(1-1/m))
						out
					}, # Close delta
					density = function(n, param, model){
						mu <- param[, 1]
						sigma <- exp(param[,2])
						exp.z.mu.sig <- exp(-(n-mu)/sigma)
						exp(-exp.z.mu.sig) * exp.z.mu.sig /sigma
					},
					rng = function(n, param, model){
						u <- runif(n)
						param[, 1] - exp(param[,2]) * log(-log(u))
					}, # end rng
					prob = function(n, param, model){
						exp( -exp(-(n-param[, 1])/exp(param[,2])))
					},
					quant = function(n, param, model){
						param[, 1] - exp(param[,2]) * log(-log(n))
					}
)
