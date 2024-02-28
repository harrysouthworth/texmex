# no covariates

rain.fit <- evm(rain,th=30)
par(mfrow=c(1,1))
texmex:::plotrl.evmOpt(rain.fit,RetPeriodRange=c(1,2000),main="Coles (2001) figure 4.5\nGPD Return Level Plot")

sealevel.fit <- evm(portpirie$SeaLevel,family=gev)
texmex:::plotrl.evmOpt(sealevel.fit,main="Coles (2001), Figure 3.5\nGEV Return Level Plot")

# with covariates

for(Family in list(gpd,gev)){
  pst <- function(msg) texmex:::texmexPst(msg,Family=Family)
  set.seed(20130513)
  n <- 100
  X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
  param <- switch(Family$name,GPD=X,GEV=cbind(5,X))
  th <- switch(Family$name,GPD=0,GEV=-Inf)
  X$Y <- Family$rng(n,param,list(threshold=th))
  fit <- evm(Y,data=X,phi=~a,xi=~b,th=th,family=Family)

  expect_error(texmex:::plotrl.evmOpt(fit), info=pst("plotrl.evmOpt: failure for model with covariates"))
}
