test.bootmex <- function(){ # this is a weak test - it tests the structure
  # of the output but not the correctness of the bootstrap coefficients; it will
  # also catch ERRORs (as opposed to FAILUREs) if the code breaks.  For strong
  # testing of this function, run test.predict.mex
  
  set.seed(20120327)
  smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
  wmarmod <- migpd(winter, mqu=.7,  penalty="none")
  
  mySdep <- mexDependence(smarmod, dqu=.7, which=1)
  myWdep <- mexDependence(wmarmod, dqu=.7, which=1)
  
  R <- 20
  
  mySboot <- bootmex(mySdep, R=R,trace=R+1)
  myWboot <- bootmex(myWdep, R=R,trace=R+1)
  
  checkEqualsNumeric(coef(mySdep)[[2]], mySboot$simpleDep, msg="bootmex: summer simpleDep from call with mar model")
  checkEqualsNumeric(coef(myWdep)[[2]], myWboot$simpleDep, msg="bootmex: winter simpleDep from call with mar model")
  
  checkEqualsNumeric(coef(smarmod), coef(mySboot$simpleMar), msg="bootmex: summer simpleMar from call with mar model")
  checkEqualsNumeric(coef(wmarmod), coef(myWboot$simpleMar), msg="bootmex: winter simpleMar from call with mar model")
  
  checkEqualsNumeric(dim(mySboot$boot[[1]]$Z)[2], dim(mySdep$dependence$Z)[2], msg="bootmex: summer dim of residuals from call with mar model")
  checkEqualsNumeric(dim(myWboot$boot[[1]]$Z)[2], dim(myWdep$dependence$Z)[2], msg="bootmex: winter dim of residuals from call with mar model")
  
  checkEqualsNumeric(dim(mySboot$boot[[1]]$dependence), dim(coef(mySdep)[[2]]), msg="bootmex: summer dim of coefficients from call with mar model")
  checkEqualsNumeric(dim(myWboot$boot[[1]]$dependence), dim(coef(myWdep)[[2]]), msg="bootmex: winter dim of coefficients from call with mar model")
  
  checkEqualsNumeric(R, length(mySboot$boot),msg="bootmex: number of bootstrap samples")
  checkEqualsNumeric(R, length(myWboot$boot),msg="bootmex: number of bootstrap samples")
  
  checkEqualsNumeric(dim(summer),dim(mySboot$boot[[1]]$Y),msg="bootmex: size of bootstrap data set")
  checkEqualsNumeric(dim(winter),dim(myWboot$boot[[5]]$Y),msg="bootmex: size of bootstrap data set")
  
  smexmod <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="gumbel",constrain=FALSE, which=1)
  wmexmod <- mex(winter, mqu=.7,  penalty="none", margins="gumbel",constrain=FALSE, which=1)
  
  mySboot <- bootmex(smexmod, R=R,trace=R+1)
  myWboot <- bootmex(wmexmod, R=R,trace=R+1)
  
  checkEqualsNumeric(coef(smexmod)[[2]], mySboot$simpleDep, msg="bootmex: summer simpleDep from call with mex model")
  checkEqualsNumeric(coef(smexmod)[[1]], coef(mySboot$simpleMar), msg="bootmex: summer simpleMar from call with mex model")
  
  checkEqualsNumeric(coef(wmexmod)[[2]], myWboot$simpleDep, msg="bootmex: winter simpleDep from call with mex model")
  checkEqualsNumeric(coef(wmexmod)[[1]], coef(myWboot$simpleMar), msg="bootmex: winter simpleMar from call with mex model")
  
  checkEqualsNumeric(R, length(mySboot$boot),msg="bootmex: number of bootstrap samples")
  checkEqualsNumeric(R, length(myWboot$boot),msg="bootmex: number of bootstrap samples")
  
  checkEqualsNumeric(dim(summer),dim(mySboot$boot[[1]]$Y),msg="bootmex: size of bootstrap data set")
  checkEqualsNumeric(dim(winter),dim(myWboot$boot[[5]]$Y),msg="bootmex: size of bootstrap data set")
  
  smexmod.1 <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="laplace",constrain=FALSE, which=1)
  smexmod.2 <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="laplace",constrain=TRUE,v=2, which=1)
  mySboot.1 <- bootmex(smexmod.1,R=R,trace=R+1)
  mySboot.2 <- bootmex(smexmod.2,R=R,trace=R+1)
  
  checkEqualsNumeric(coef(smexmod.1)[[2]], mySboot.1$simpleDep, msg="bootmex: summer simpleDep from call with mex model, constrain=FASLE")
  checkEqualsNumeric(coef(smexmod.1)[[1]], coef(mySboot.1$simpleMar), msg="bootmex: summer simpleMar from call with mex model, constrain=FASLE")
  checkEqualsNumeric(coef(smexmod.2)[[2]], mySboot.2$simpleDep, msg="bootmex: summer simpleDep from call with mex model, v=2")
  checkEqualsNumeric(coef(smexmod.2)[[1]], coef(mySboot.2$simpleMar), msg="bootmex: summer simpleMar from call with mex model, v=2")
  
  # check execution of for 2-d data
  
  wavesurge.fit <- migpd(wavesurge,mqu=.7)
  wavesurge.mex <- mexDependence(wavesurge.fit, dqu=0.8,which=1)
  R <- 20
  
  wavesurge.boot <- bootmex(wavesurge.mex,R=R,trace=R+1)
  
  checkEqualsNumeric(dim(wavesurge.boot$boot[[1]]$Z)[2],1,msg="bootmex: execution for 2-d data")
  checkEquals(dimnames(wavesurge.boot$boot[[1]]$Z)[[2]],names(wavesurge)[2],msg="bootmex: execution for 2-d data")
  checkEqualsNumeric(length(wavesurge.boot$boot),R,msg="bootmex: execution for 2-d data")
}

test.chi <- function(){
  
  # independent implementation of chi and chibar, Janet Heffernan personal code library
  
  
  .ChiFunction <- function(data, nLevels){
    .Cfunction <- function(data, nLevels){
      rowWiseMax <- apply(data, 1, max)
      rowWiseMin <- apply(data, 1, min)
      u <- seq(min(rowWiseMax) + 1/(2 * nLevels), max(rowWiseMin) - 1/(2 * nLevels), length = nLevels)
      Cu <- sapply(1:nLevels,function(i, rmax, u){ mean(rmax < u[i]) }, rmax = rowWiseMax, u=u)
      CbarU <- sapply(1:nLevels,function(i, rmin, u){ mean(rmin > u[i]) }, rmin=rowWiseMin, u=u)
      list(u = u, Cu = Cu, CbarU = CbarU)
    } 
    TransUniform <- function(x){
      .transUniform <- function(x){
        if (is.R()){
          rank(x,ties.method="first") / (length(x) + 1) # original version
        }
        else {
          rank(x) / (length(x) + 1) # original version
        }
      }
      
      
      if(length(dim(x)) > 0)apply(x,2,.transUniform)
      else .transUniform(x)
    } 
    C <- .Cfunction(TransUniform(data), nLevels = nLevels)
    u <- C$u
    Cu <- C$Cu
    CbarU <- C$CbarU
    ChiU <- 2 - log(Cu)/log(u)
    ChiBarU <- (2 * log(1 - u))/log(CbarU) - 1
    n <- nrow(data)
    
    #variances of chi and chibar
    varChi <- ((1/log(u)^2 * 1)/Cu^2 * Cu * (1 - Cu))/n
    varChiBar <- (((4 * log(1 - u)^2)/(log(CbarU)^4 * CbarU^2)) * CbarU * (1 - CbarU))/n
    
    #upper and lower 95% conf int bounds for chi and chibar; these are based on normal approx with further functional constraints imposed
    z.975 <- qnorm(1 - 0.05/2)
    ChiLower <- ChiU - z.975 * sqrt(varChi)
    ChiUpper <- ChiU + z.975 * sqrt(varChi)
    
    ChiLbound <- numeric(length(u))
    ChiLbound[u>0.5] <- 2 - log(2 *u[u > 0.5] - 1)/log(u[u > 0.5])
    ChiLbound[u<=0.5] <- -Inf
    
    ChiLower <- apply(cbind(ChiLower, ChiLbound), 1, max)
    ChiUpper[ChiUpper > 1] <- 1
    
    ChiBarLower <- ChiBarU - z.975 * sqrt(varChiBar)
    ChiBarUpper <- ChiBarU + z.975 * sqrt(varChiBar)
    ChiBarLower[ChiBarLower < -1] <- -1
    ChiBarUpper[ChiBarUpper > 1] <- 1
    
    list(u = C$u, Cu = C$Cu, CbarU = C$CbarU, 
         Chi = ChiU, ChiBar = ChiBarU, 
         ChiLower = ChiLower, ChiUpper = ChiUpper,
         ChiBarLower = ChiBarLower, ChiBarUpper = ChiBarUpper,
         n = n)
  }
  
  
  
  #*************************************************************
  
  nq <- 1000
  chi.JH <- .ChiFunction(wavesurge,nLevels=nq)
  chi <- chi(wavesurge,nq=nq,qlim=range(chi.JH$u),trunc= TRUE)
  
  checkEqualsNumeric(chi.JH$u,chi$quantile,msg="chi: u")
  checkEqualsNumeric(chi.JH$Chi,chi$chi[,2],msg="chi: Chi")
  checkEqualsNumeric(chi.JH$ChiLower,chi$chi[,1],msg="chi: ChiLower")
  checkEqualsNumeric(chi.JH$ChiUpper,chi$chi[,3],msg="chi: ChiUpper")
  checkEqualsNumeric(chi.JH$ChiBar, chi$chibar[,2],msg="chi: ChiBar")
  checkEqualsNumeric(chi.JH$ChiBarLower, chi$chibar[,1],msg="chi: ChiBarLower")
  checkEqualsNumeric(chi.JH$ChiBarUpper, chi$chibar[,3],msg="chi: ChiBarUpper")
}

test.plot.chi <- function(){
  chi <- chi(wavesurge)
  par(mfrow=c(1,2),pty="m")
  res <- plot(chi,mainChiBar="Figure 8.11 of Coles (2001)\nChi Bar")
  plot(chi, show=c("Chi"=FALSE,"ChiBar"=TRUE))
  plot(chi, show=c("Chi"=TRUE,"ChiBar"=FALSE))
  checkEquals(res,NULL,msg = "plot.chi: check successful execution")
}

test.copula <- function(){
  fun <- function(d) apply(d,2,function(x)(1:n)[rank(x)])/(1+n)
  n <- 200
  
  u2 <- cbind(sample(n),sample(n))
  d2 <- fun(u2)
  
  u3 <- cbind(sample(n),sample(n),sample(n))
  d3 <- fun(u3)
  
  checkEqualsNumeric(d2,copula(u2)$copula,msg="copula: 2 dimensional")
  checkEqualsNumeric(d3,copula(u3)$copula,msg="copula: 3 dimensional")
  
  op <- options()
  options(show.error.messages=FALSE)
  checkException(copula(TRUE),msg="copula: exception")
  checkException(copula("text"),msg="copula: exception")
  options(op)
}

test.dgpd <- function(){
  evd.dgpd <- texmex:::.evd.dgpd
  
  myTest <- function(sig,xi,thresh,msg){
    myd <- sapply(1:nreps,function(i) dgpd(x[,i], sig[i], xi[i],u=thresh[i]))
    ed <- sapply(1:nreps, function(i) evd.dgpd(x[,i], loc=thresh[i], scale=sig[i], shape=xi[i]))
    checkEqualsNumeric(ed,myd,msg=msg)
  }
  
  set.seed(20101111)
  
  #*************************************************************
  # 6.12. Test dgpd. Note that .evd.dgpd is NOT vectorized.
  
  nreps <- 100
  nsim <- 1000
  p <- matrix(runif(2*nreps, -1, 1),ncol=2)
  p[, 1] <- p[, 1] + 1
  thresh <- rep(0,nreps)
  
  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))
  
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="dgpd: random xi")
  
  #*************************************************************
  # 6.13. Test dgpd when some or all of xi == 0
  
  p[sample(1:nreps,nreps/2),2] <- 0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="dgpd: some zero xi")
  
  p[,2] <-  0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="dgpd: all zero xi")
  
  #*************************************************************
  # 6.14. Test vectorization of dgpd.
  
  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)
  thresh <- rnorm(nsim)
  
  x <- rgpd(nsim, sig, xi,u=thresh)
  myd <- dgpd(x, sig, xi,u=thresh)
  
  ed <- sapply(1:nsim, function(i) evd.dgpd(x[i], loc=thresh[i], scale=sig[i], shape=xi[i]))
  checkEqualsNumeric(ed,myd,msg="dgpd: vectorisation")
  
  #*************************************************************
  # 6.15 test log.d argument
  
  ld <- dgpd(x,sig,xi,u=thresh,log.d=TRUE)
  checkEqualsNumeric(myd,exp(ld),msg="dgpd: log density")
}

test.closures <- function() {
  make.mvn.prior <- texmex:::.make.mvn.prior
  make.quad.prior <- texmex:::.make.quadratic.penalty
  make.spd.matrix <- texmex:::.random.spd.matrix
  for (count in 1:100) {
    dimn <- sample(10, size=1)
    cov.matrix <- make.spd.matrix(dimn)
    centre <- rexp(dimn)
    mvnprior <- make.mvn.prior(list(centre, cov.matrix))
    point <- rexp(dimn)
    checkEqualsNumeric(mvnprior(point),
                       dmvnorm(point, centre, cov.matrix, log=TRUE),
                       "efficient.closures: multivariate Gaussian prior")
    quadprior <- make.quad.prior(list(centre, cov.matrix))
    checkEqualsNumeric(quadprior(point),
                       mahalanobis(point, centre, cov.matrix),
                       "efficient.closures: Mahalanobis distance")
    # or "A-norm", as it's otherwise called
  }
}

test.endPoint <- function(){
  
  for(Family in list(gpd,gev)){
    pst <- function(msg) texmexPst(msg,Family=Family)
    set.seed(20130617)
    
    for(i in 1:5){
      th <- switch(Family$name,GPD=0.3,GEV=-Inf)
      fit <- evm(rnorm(100),th=th,family=Family)
      co <- coef(fit)
      ep.current <- endPoint(fit,verbose=FALSE,.unique=TRUE)
      ep.target <- switch(Family$name,GPD=ifelse(co[2] < 0, th-exp(co[1])/co[2],Inf),
                          GEV=ifelse(co[3] < 0, co[1]-exp(co[2])/co[3],Inf))
      checkEqualsNumeric(ep.current, ep.target, msg=pst("endPoint: check calc for evmOpt no covariates"))
      
      fit <- evm(rnorm(100),th=th,family=Family,method="simulate",trace=50000)
      co <- coef(fit$map)
      ep.current <- endPoint(fit,verbose=FALSE,.unique=TRUE)
      ep.target <- switch(Family$name,GPD=ifelse(co[2] < 0, th-exp(co[1])/co[2],Inf), 
                          GEV=ifelse(co[3] < 0, co[1]-exp(co[2])/co[3],Inf))
      checkEqualsNumeric(ep.current, ep.target, msg=pst("endPoint: check calc for evmSim no covariates"))
    }  
    
    # test with covariates
    
    n <- 50
    mu <- 1
    for(i in 1:5){
      X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
      param <- switch(Family$name,GPD=X,GEV=cbind(mu,X))
      th <- switch(Family$name,GPD=0,GEV=-Inf)
      X$Y <- Family$rng(n,param,list(threshold=th))
      fit <- evm(Y,data=X,phi=~a,xi=~b,th=th,family=Family)
      lp <- linearPredictors(fit)$link
      ep.current <- endPoint(fit,verbose=FALSE,.unique=FALSE)
      ep.target <- switch(Family$name,GPD=ifelse(lp[,2] < 0, th-exp(lp[,1])/lp[,2],Inf),
                          GEV=ifelse(lp[,3] < 0, lp[,1]-exp(lp[,2])/lp[,3],Inf))
      checkEqualsNumeric(ep.current, ep.target, msg=pst("endPoint: check calc for evmSim with covariates"))
    }
  }
}

test.evm <- function(){
  tol <- 0.01
  
  ###################################################################
  # 1.3 Reproduce loglik, parameter estimates and covariance on page 85
  #    of Coles. Will not be exact because fitting routines differ:
  #    gpd works with phi=log(sigma). Can only reproduce cell [2,2]
  #    of covariance.
  
  cparas <- c(7.44, 0.184)
  cse <- c(0.958432, 0.101151)
  
  ccov <- matrix(c(.9188, -.0655, -.0655, .0102), nrow=2)
  cloglik <- -485.1
  
  mod <- evm(rain, th=30, penalty="none")
  mod.coef <- coef(mod)
  
  mod.coef[1] <- exp(mod.coef[1])
  names(mod.coef)[1] <- "sigma"
  
  mod.loglik <- mod$loglik
  mod.cov22 <- mod$cov[2, 2]
  
  checkEqualsNumeric(cparas, mod.coef, tolerance = tol,
                     msg="gpd: parameter ests page 85 Coles")
  checkEqualsNumeric(cse[2], mod$se[2], tolerance = tol,
                     msg="gpd: standard errors page 85 Coles")
  checkEqualsNumeric(ccov[2,2], mod$cov[2, 2], tolerance = tol,
                     msg="gpd: Covariance page 85 Coles")
  checkEqualsNumeric(cloglik, mod.loglik, tolerance = tol,
                     msg="gpd: loglik page 85 Coles")
  
  ###################################################################
  #   Logical checks on the effect of Gaussian penalization. The smaller the
  #    variance, the more the parameter should be drawn towards the
  #    mean.
  
  # 2.1 Tests for xi being drawn to 0
  
  gp1 <- list(c(0, 0), diag(c(10^4, .25)))
  gp2 <- list(c(0, 0), diag(c(10^4, .05)))
  
  mod1 <- evm(rain, th=30, priorParameters=gp1)
  mod2 <- evm(rain, th=30, priorParameters=gp2)
  
  checkTrue(coef(mod)[2] > coef(mod1)[2],
            msg="gpd: Gaussian penalization xi being drawn to 0")
  checkTrue(coef(mod1)[2] > coef(mod2)[2],
            msg="gpd: Gaussian penalization xi being drawn to 0")
  
  # 2.2 Tests for phi being drawn to 0
  
  gp3 <- list(c(0, 0), diag(c(1, 10^4)))
  gp4 <- list(c(0, 0), diag(c(.1, 10^4)))
  
  mod3 <- evm(rain, th=30, priorParameters=gp3)
  mod4 <- evm(rain, th=30, priorParameters=gp4)
  
  checkTrue(coef(mod)[1] > coef(mod3)[1],
            msg="gpd: Gaussian penalization phi being drawn to 0")
  checkTrue(coef(mod3)[1] > coef(mod4)[1],
            msg="gpd: Gaussian penalization phi being drawn to 0")
  
  # 2.3 Tests for xi being drawn to 1
  gp5 <- list(c(0, 1), diag(c(10^4, .25)))
  gp6 <- list(c(0, 1), diag(c(10^4, .05)))
  
  mod5 <- evm(rain, th=30, priorParameters=gp5)
  mod6 <- evm(rain, th=30, priorParameters=gp6)
  
  checkTrue(1 - coef(mod)[2] > 1 - coef(mod5)[2],
            msg="gpd: Gaussian penalization xi being drawn to 1")
  checkTrue(1 - coef(mod1)[2] > 1 - coef(mod6)[2],
            msg="gpd: Gaussian penalization xi being drawn to 1")
  
  # 2.4 Tests for phi being drawn to 4 (greater than mle for phi)
  
  gp7 <- list(c(4, 0), diag(c(1, 10^4)))
  gp8 <- list(c(4, 0), diag(c(.1, 10^4)))
  
  mod7 <- evm(rain, th=30, priorParameters=gp7)
  mod8 <- evm(rain, th=30, priorParameters=gp8)
  
  checkTrue(4 - coef(mod)[1] > 4 - coef(mod7)[1],
            msg="gpd: Gaussian penalization phi being drawn to 4")
  checkTrue(4 - coef(mod3)[1] > 4 - coef(mod8)[1],
            msg="gpd: Gaussian penalization phi being drawn to 4")
  
  ###################################################################
  #   Logical checks on the effect of penalization using lasso or L1 penalization. The smaller the
  #    variance, the more the parameter should be drawn towards the
  #    mean.
  
  # 2a.1 Tests for xi being drawn to 0
  
  gp1 <- list(c(0, 0), solve(diag(c(10^4, .25))))
  gp2 <- list(c(0, 0), solve(diag(c(10^4, .05))))
  
  mod1 <- evm(rain, th=30, priorParameters=gp1, penalty="lasso")
  mod2 <- evm(rain, th=30, priorParameters=gp2, penalty="lasso")
  
  checkTrue(coef(mod)[2] > coef(mod1)[2],
            msg="gpd: lasso penalization xi being drawn to 0")
  checkTrue(coef(mod1)[2] > coef(mod2)[2],
            msg="gpd: lasso penalization xi being drawn to 0")
  
  # 2a.2 Tests for phi being drawn to 0
  
  gp3 <- list(c(0, 0), solve(diag(c(1, 10^4))))
  gp4 <- list(c(0, 0), solve(diag(c(.1, 10^4))))
  
  mod3 <- evm(rain, th=30, priorParameters=gp3, penalty="lasso")
  mod4 <- evm(rain, th=30, priorParameters=gp4, penalty="lasso")
  
  checkTrue(coef(mod)[1] > coef(mod3)[1],
            msg="gpd: lasso penalization phi being drawn to 0")
  checkTrue(coef(mod3)[1] > coef(mod4)[1],
            msg="gpd: lasso penalization phi being drawn to 0")
  
  # 2a.3 Tests for xi being drawn to 1
  gp5 <- list(c(0, 1), solve(diag(c(10^4, .25))))
  gp6 <- list(c(0, 1), solve(diag(c(10^4, .05))))
  
  mod5 <- evm(rain, th=30, priorParameters=gp5, penalty="lasso")
  mod6 <- evm(rain, th=30, priorParameters=gp6, penalty="lasso")
  
  checkTrue(1 - coef(mod)[2] > 1 - coef(mod5)[2],
            msg="gpd: lasso penalization xi being drawn to 1")
  checkTrue(1 - coef(mod1)[2] > 1 - coef(mod6)[2],
            msg="gpd: lasso penalization xi being drawn to 1")
  
  # 2a.4 Tests for phi being drawn to 4 (greater than mle for phi)
  
  gp7 <- list(c(4, 0), solve(diag(c(1, 10^4))))
  gp8 <- list(c(4, 0), solve(diag(c(.1, 10^4))))
  
  mod7 <- evm(rain, th=30, priorParameters=gp7, penalty="lasso")
  mod8 <- evm(rain, th=30, priorParameters=gp8, penalty="lasso")
  
  checkTrue(4 - coef(mod)[1] > 4 - coef(mod7)[1],
            msg="gpd: lasso penalization phi being drawn to 4")
  checkTrue(4 - coef(mod3)[1] > 4 - coef(mod8)[1],
            msg="gpd: lasso penalization phi being drawn to 4")
  
  ########################################################
  # Tests on including covariates. Once more, gpd.fit in ismev
  # works with sigma inside the optimizer, so we need to tolerate
  # some differences and standard errors might be a bit out.
  
  # 3.0 Reproduce Coles, page 119. Reported log-likelihood is -484.6.
  
  rtime <- (1:length(rain))/1000
  d <- data.frame(rainfall = rain, time=rtime)
  
  mod <- evm(rainfall, th=30, data=d, phi= ~ time, penalty="none")
  
  checkEqualsNumeric(-484.6, mod$loglik, tolerance = tol,
                     msg="gpd: loglik Coles page 119")
  
  ####################################################################
  # 3.1 Use liver data, compare with ismev.
  #     These are not necessarily sensible models!
  #     Start with phi alone.
  
  mod <- evm(ALT.M, qu=.7, data=liver,
             phi = ~ ALT.B + dose, xi = ~1,
             penalty="none", cov="observed")
  
  m <- model.matrix(~ ALT.B + dose, liver)
  
  ismod <- texmex:::.ismev.gpd.fit(liver$ALT.M,
                                   threshold=quantile(liver$ALT.M, .7),
                                   ydat = m, sigl=2:ncol(m),
                                   siglink=exp, show=FALSE)
  
  checkEqualsNumeric(ismod$mle,
                     coef(mod),
                     tolerance=tol,
                     msg="gpd: covariates in phi only, point ests")
  
  # SEs for phi will not be same as for sigma, but we can test xi
  checkEqualsNumeric(ismod$se[length(ismod$se)],
                     mod$se[length(mod$se)],
                     tolerance=tol,
                     msg="gpd: covariates in phi only, standard errors")
  
  ######################################################################
  # 3.2 Test xi alone.
  mod <- evm(log(ALT.M / ALT.B), qu=.7, data=liver,
             phi = ~1, xi = ~ ALT.B + dose,
             penalty="none")
  
  m <- model.matrix(~ ALT.B + dose, liver)
  
  ismod <- texmex:::.ismev.gpd.fit(log(liver$ALT.M / liver$ALT.B),
                                   threshold=quantile(log(liver$ALT.M / liver$ALT.B), .7),
                                   ydat = m, shl=2:ncol(m), show=FALSE)
  mco <- coef(mod)
  mco[1] <- exp(mco[1])
  
  checkEqualsNumeric(ismod$mle,
                     mco,
                     tolerance=tol,
                     msg="gpd: covariates in xi only: point ests")
  # SEs for phi will not be same as for sigma, but we can test xi
  checkEqualsNumeric(ismod$se[-1],
                     mod$se[-1],
                     tolerance = tol,
                     msg="gpd: covariates in xi only: standard errors")
  
  ######################################################################
  # 3.3 Test phi & xi simultaneously. Use simulated data.
  
  set.seed(25111970)
  
  makeDataGpd <- function(a,b,n=500,u=10)
    # lengths of a and b should divide n exactly
    # returns data set size 2n made up of uniform variates (size n) below threshold u and
    # gpd (size n) with scale parameter exp(a) and shape b above threshold u
  {
    gpd <- rgpd(n,exp(a),b,u=u)
    unif <- runif(n,u-10,u)
    as.data.frame(cbind(a=a,b=b,y=c(gpd,unif)))
  }
  
  mya <- seq(0.1,1,len=10)
  myb <- rep(c(-0.2,0.2),each=5)
  data <- makeDataGpd(mya,myb)
  m <- model.matrix(~ a+b, data)
  
  mod <- evm(y,qu=0.7,data=data,phi=~a,xi=~b,penalty="none")
  ismod <- texmex:::.ismev.gpd.fit(data$y,
                                   threshold=quantile(data$y,0.7),
                                   ydat=m,shl=3,sigl=2,
                                   siglink=exp,
                                   show=FALSE)
  
  checkEqualsNumeric(ismod$mle,
                     coef(mod),
                     tolerance = tol,
                     msg="gpd: covariates in phi and xi: point ests")
  checkEqualsNumeric(ismod$se,
                     sqrt(diag(mod$cov)),
                     tolerance = tol,
                     msg="gpd: covariates in phi and xi: std errs")
  
  ####################################################################
  # Check that using priors gives expected behaviour when covariates are included.
  
  # 2.1 Tests for xi being drawn to 0
  
  myb <- rep(c(0.5,1.5),each=5)
  data <- makeDataGpd(a=1,b=myb,n=3000)
  
  gp1 <- list(c(0, 0, 0), diag(c(10^4, 0.25, 0.25)))
  gp2 <- list(c(0, 0, 0), diag(c(10^4, 0.25, 0.01)))
  
  mod0 <- evm(y,qu=0.6,data=data,phi=~1,xi=~b,penalty="none")
  mod1 <- evm(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp1)
  mod2 <- evm(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp2)
  
  checkTrue(all(abs(coef(mod0)[2:3]) > abs(coef(mod1)[2:3])),
            msg="gpd: with covariates, xi drawn to zero")
  checkTrue(abs(coef(mod1)[3]) > abs(coef(mod2)[3]),
            msg="gpd: with covariates, xi drawn to zero")
  
  # 2.2 Tests for phi being drawn to 0
  
  # HS. Changed a to mya due to scoping problems in S+. The issue is very general
  # and affects (for example) lm(~a, data, method="model.frame"), so it's kind of
  # by design.
  mya <- seq(0.1,1,len=10)
  data <- makeDataGpd(-3 + mya,b=-0.1,n=3000)
  data$a <- rep(mya, len=nrow(data))
  
  gp4 <- list(c(0, 0, 0), diag(c(1, 1, 10^4)))
  gp5 <- list(c(0, 0, 0), diag(c(0.1, 0.1, 10^4)))
  
  mod3 <- evm(y,qu=0.6,data=data,phi=~a,xi=~1,penalty="none")
  mod4 <- evm(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp4)
  mod5 <- evm(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp5)
  
  checkTrue(all(abs(coef(mod3)[1:2]) > abs(coef(mod4)[1:2])),
            msg="gpd: with covariates, phi being drawn to 0")
  checkTrue(all(abs(coef(mod4)[1:2]) > abs(coef(mod5)[1:2])),
            msg="gpd: with covariates, phi being drawn to 0")
  
  # 2.3 Tests for xi being drawn to 2
  myb <- rep(c(-0.5,0.5),each=5)
  data <- makeDataGpd(a=1,b=myb,n=3000)
  
  gp7 <- list(c(0, 2, 2), diag(c(10^4, 0.25, 0.25)))
  gp8 <- list(c(0, 2, 2), diag(c(10^4, 0.05, 0.05)))
  
  mod6 <- evm(y,qu=0.6,data=data,phi=~1,xi=~b,penalty="none")
  mod7 <- evm(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp7)
  mod8 <- evm(y,qu=0.6,data=data,phi=~1,xi=~b,priorParameters=gp8)
  
  checkTrue(all(abs(2 - coef(mod6)[2:3]) > abs(2 - coef(mod7)[2:3])),
            msg="gpd: with covariates, xi drawn to 2")
  checkTrue(all(abs(2 - coef(mod7)[2:3]) > abs(2 - coef(mod8)[2:3])),
            msg="gpd: with covariates, xi drawn to 2")
  
  # 2.4 Tests for phi being drawn to 4
  
  mya <- seq(0.1,1,len=10)
  data <- makeDataGpd(2 + mya,b=-0.1,n=3000)
  data$a <- rep(mya, len=nrow(data))
  
  gp10 <- list(c(0, 4, 0), diag(c(10^4, 1,   10^4)))
  gp11 <- list(c(0, 4, 0), diag(c(10^4, 0.1, 10^4)))
  
  mod9 <- evm(y,qu=0.6,data=data,phi=~a,xi=~1,penalty="none")
  mod10 <- evm(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp10)
  mod11 <- evm(y,qu=0.6,data=data,phi=~a,xi=~1,priorParameters=gp11)
  
  checkTrue(abs(4 - coef(mod9)[2])  > abs(4 - coef(mod10)[2]),
            msg="gpd: with covariates, phi drawn to 4")
  checkTrue(abs(4 - coef(mod10)[2]) > abs(4 - coef(mod11)[2]),
            msg="gpd: with covariates, phi drawn to 4")
  
  #*************************************************************
  # 4.1. Test reproducibility
  set.seed(20101110)
  save.seed <- .Random.seed
  
  set.seed(save.seed)
  bmod <- evm(ALT.M, data=liver,
              th=quantile(liver$ALT.M, .7),
              iter=1000, thin=1, verbose=FALSE, method="sim")
  
  set.seed(save.seed)
  bmod2 <- evm(ALT.M, data=liver,
               th=quantile(liver$ALT.M, .7),
               iter=1000, thin=1, verbose=FALSE, method="sim")
  
  checkEqualsNumeric(bmod$param, bmod2$param,
                     msg="gpd: test simulation reproducibility 1")
  
  set.seed(bmod$seed)
  bmod3 <- evm(ALT.M, data=liver,
               th=quantile(liver$ALT.M, .7),
               iter=1000, thin=1, verbose=FALSE, method="sim")
  checkEqualsNumeric(bmod$param, bmod3$param,
                     msg="gpd: test simulation reproducibility 2")
  
  #*************************************************************
  # 4.2. Logical test of burn-in
  
  checkEqualsNumeric(nrow(bmod$chains) - bmod$burn, nrow(bmod$param),
                     msg="gpd: Logical test of burn-in 1")
  
  iter <- sample(500:1000,1)
  burn <- sample(50,1)
  bmod2 <- evm(ALT.M, data=liver, th=quantile(liver$ALT.M, .7),
               iter=iter, burn=burn, thin=1, verbose=FALSE, method="sim")
  
  checkEqualsNumeric(iter-burn, nrow(bmod2$param),
                     msg="gpd: Logical test of burn-in 2")
  
  #*************************************************************
  # 4.3. Logical test of thinning
  
  thin <- 0.5
  iter <- 1000
  bmod <- evm(ALT.M, data=liver, th=quantile(liver$ALT.M, .7),
              iter=iter, thin = thin,verbose=FALSE, method="sim")
  
  checkEqualsNumeric((nrow(bmod$chains) - bmod$burn) * thin, nrow(bmod$param),
                     msg="gpd: Logical test of thinning 1")
  
  thin <- 2
  iter <- 1000
  bmod <- evm(ALT.M, data=liver, th=quantile(liver$ALT.M, .7),
              iter=iter, thin = thin, verbose=FALSE, method="sim")
  
  checkEqualsNumeric((nrow(bmod$chains) - bmod$burn) / thin, nrow(bmod$param),
                     msg="gpd: Logical test of thinning 1")
  
  #*************************************************************
  # 5.1. Test of gev family: point estimates and cov matrices
  coles <- c(3.87, .198, -.050) # From page 59 of Coles
  mOpt <- evm(SeaLevel, data=portpirie, family=gev, penalty="none")
  mSim <- evm(SeaLevel, data=portpirie, family=gev, method="sim",trace=100000)
  coOpt <- coef(mOpt)
  coSim <- coef(mSim)
  coOpt[2] <- exp(coOpt[2])
  coSim[2] <- exp(coSim[2])
  # check point estimates
  checkEqualsNumeric(coles, coOpt, tolerance = tol,
                     msg="gev: optimisation, parameter ests page 59 Coles")
  checkEqualsNumeric(coles, coSim, tolerance = tol,
                     msg="gev: simulation, parameter ests page 59 Coles")
  
  # Check non-sigma elements of covariance
  coles <- matrix(c(.000780, -.00107,
                    -.00107, .00965), ncol=2)
  coOpt <- mOpt$cov[c(1,3), c(1,3)]
  coSim <- cov(mSim$param)[c(1,3), c(1,3)]
  checkEqualsNumeric(coles, coOpt, tolerance=tol,
                     msg="gev: optimisation, covariance page 59 coles")
  checkEqualsNumeric(coles, coSim, tolerance=tol,
                     msg="gev: simulation, covariance page 59 coles")
  mcOpt <- max(abs(coles - coOpt))
  mcSim <- max(abs(coles - coSim))
  checkEqualsNumeric(0, mcOpt, tolerance=tol, msg="gev: optimisation, max abs covariance")
  checkEqualsNumeric(0, mcSim, tolerance=tol, msg="gev: simulation, max abs covariance")
  
  # Check log-likelihood
  coles <- 4.34
  co <- mOpt$loglik
  checkEqualsNumeric(coles, co, tolerance=tol,
                     msg="gev: optimisation, loglik page 59 coles")
  
  ###################################################################
  # 5.2   GEV - Logical checks on the effect of Gaussian penalization. The smaller the
  #    variance, the more the parameter should be drawn towards the
  #    mean.
  
  # Tests for xi being drawn to 0
  
  gp1 <- list(c(0, 0, 0), diag(c(10^4, 10^4, .25)))
  gp2 <- list(c(0, 0, 0), diag(c(10^4, 10^4, .05)))
  
  mod1 <- evm(SeaLevel, data=portpirie, family=gev, priorParameters=gp1)
  mod2 <- evm(SeaLevel, data=portpirie, family=gev, priorParameters=gp2)
  
  checkTrue(abs(coef(mOpt)[3]) > abs(coef(mod1)[3]),
            msg="gev: Gaussian penalization xi being drawn to 0")
  checkTrue(abs(coef(mOpt)[3]) > abs(coef(mod2)[3]),
            msg="gev: Gaussian penalization xi being drawn to 0")
  checkTrue(abs(coef(mod1)[3]) > abs(coef(mod2)[3]),
            msg="gev: Gaussian penalization xi being drawn to 0")
  
  # Tests for phi being drawn to 0
  
  gp3 <- list(c(0, 0, 0), diag(c(10^4, 1, 10^4)))
  gp4 <- list(c(0, 0, 0), diag(c(10^4, .1, 10^4)))
  
  mod3 <- evm(SeaLevel, data=portpirie, family=gev, priorParameters=gp3)
  mod4 <- evm(SeaLevel, data=portpirie, family=gev, priorParameters=gp4)
  
  checkTrue(abs(coef(mOpt)[2]) > abs(coef(mod3)[2]),
            msg="gev: Gaussian penalization phi being drawn to 0")
  checkTrue(abs(coef(mOpt)[2]) > abs(coef(mod4)[2]),
            msg="gev: Gaussian penalization phi being drawn to 0")
  checkTrue(abs(coef(mod3)[2]) > abs(coef(mod4)[2]),
            msg="gev: Gaussian penalization phi being drawn to 0")
  
  # Tests for xi being drawn to 1
  gp5 <- list(c(0, 0, 1), diag(c(10^4, 10^4, .25)))
  gp6 <- list(c(0, 0, 1), diag(c(10^4, 10^4, .05)))
  
  mod5 <- evm(SeaLevel, data=portpirie, family=gev, priorParameters=gp5)
  mod6 <- evm(SeaLevel, data=portpirie, family=gev, priorParameters=gp6)
  
  checkTrue(abs(1 - coef(mOpt)[3]) > abs(1 - coef(mod5)[3]),
            msg="gev: Gaussian penalization xi being drawn to 1")
  checkTrue(abs(1 - coef(mOpt)[3]) > abs(1 - coef(mod6)[3]),
            msg="gev: Gaussian penalization xi being drawn to 1")
  checkTrue(abs(1 - coef(mod5)[3]) > abs(1 - coef(mod6)[3]),
            msg="gev: Gaussian penalization xi being drawn to 1")
  
  # Tests for phi being drawn to 0.5 (greater than mle for phi)
  
  gp7 <- list(c(0, 0.5, 0), diag(c(10^4, 1, 10^4)))
  gp8 <- list(c(0, 0.5, 0), diag(c(10^4, .1, 10^4)))
  
  mod7 <- evm(SeaLevel, data=portpirie, family=gev, priorParameters=gp7)
  mod8 <- evm(SeaLevel, data=portpirie, family=gev, priorParameters=gp8)
  
  checkTrue(abs(.5 - coef(mOpt)[2]) > abs(.5 - coef(mod7)[2]),
            msg="gpd: Gaussian penalization phi being drawn to .5")
  checkTrue(abs(.5 - coef(mOpt)[2]) > abs(.5 - coef(mod8)[2]),
            msg="gpd: Gaussian penalization phi being drawn to .5")
  
  # test above but with lasso penalty
  
  mod9 <- evm(SeaLevel, data=portpirie, family=gev, priorParameters=gp7,penalty="lasso")
  mod10 <- evm(SeaLevel, data=portpirie, family=gev, priorParameters=gp8,penalty="lasso")
  
  checkTrue(abs(.5 - coef(mOpt)[2]) > abs(.5 - coef(mod9)[2]),
            msg="gpd: Gaussian lasso penalization phi being drawn to .5")
  checkTrue(abs(.5 - coef(mOpt)[2]) > abs(.5 - coef(mod10)[2]),
            msg="gpd: Gaussian lasso penalization phi being drawn to .5")
  
  ###################################################################
  # 5.3  GEV - covariates in mu, phi and xi separately, use gev.fit from ismev package
  
  set.seed(20130614)
  makeDataGev <- function(u=0,a,b,n=500)
    # lengths of a and b should divide n exactly
    # returns data set size n distributed as GEV variates location parameter u,
    # scale parameter exp(a) and shape b
  {
    gev <- rgev(n,mu=u,exp(a),b)
    as.data.frame(cbind(u=u,a=a,b=b,y=gev))
  }
  
  myu <- rnorm(10)
  mya <- seq(0.1,1,len=10)
  myb <- rep(c(-0.2,0.2),each=5)
  
  data <- list(makeDataGev(myu,2,-0.1),
               makeDataGev(3,mya,-0.1),
               makeDataGev(3,2,myb))
  
  start1 <- c(0,1,2,-0.1)
  g1.fit <- texmex:::.ismev.gev.fit(xdat=data[[1]]$y,ydat=data[[1]],mul=1, siglink=exp,show=FALSE,muinit=start1[1:2])
  g2.fit <- texmex:::.ismev.gev.fit(xdat=data[[2]]$y,ydat=data[[2]],sigl=2,siglink=exp,show=FALSE)
  g3.fit <- texmex:::.ismev.gev.fit(xdat=data[[3]]$y,ydat=data[[3]],shl=3, siglink=exp,show=FALSE)
  t1.fit <- evm(y,mu=~u, family=gev,data=data[[1]],start=start1)
  t2.fit <- evm(y,phi=~a,family=gev,data=data[[2]])
  t3.fit <- evm(y,xi=~b, family=gev,data=data[[3]])
  
  checkEqualsNumeric(g1.fit$mle,coef(t1.fit),tolerance=tol,msg="gev: covariates in mu point est")
  checkEqualsNumeric(g2.fit$mle,coef(t2.fit),tolerance=tol,msg="gev: covariates in phi point est")
  checkEqualsNumeric(g3.fit$mle,coef(t3.fit),tolerance=tol,msg="gev: covariates in xi point est")
  
  checkEqualsNumeric(g1.fit$nllh,-t1.fit$loglik,tolerance=tol,msg="gev: covariates in mu, log-lik")
  checkEqualsNumeric(g2.fit$nllh,-t2.fit$loglik,tolerance=tol,msg="gev: covariates in phi, log-lik")
  checkEqualsNumeric(g3.fit$nllh,-t3.fit$loglik,tolerance=tol,msg="gev: covariates in xi, log-lik")
  
  checkEqualsNumeric(g1.fit$cov,t1.fit$cov,tolerance=tol, msg="gev: covariates in mu, cov")
  checkEqualsNumeric(g2.fit$cov,t2.fit$cov,tolerance=tol, msg="gev: covariates in phi, cov")
  checkEqualsNumeric(g3.fit$cov,t3.fit$cov,tolerance=tol, msg="gev: covariates in xi, cov")
  
  ######################################################################
  # 5.4 GEV - Test mu phi & xi simultaneously. Use simulated data.
  
  set.seed(25111970)
  
  data <- makeDataGev(5+2*myu,2-0.5*mya,0.1+2*myb)
  data$u <- myu
  data$a <- mya
  data$b <- myb
  
  mod <- evm(y,data=data,mu=~u,phi=~a,xi=~b,penalty="none",family=gev,start=c(5,2,2,-0.1,0,2))
  ismod <- texmex:::.ismev.gev.fit(data$y,
                                   ydat=data,mul=1,shl=3,sigl=2,
                                   siglink=exp,
                                   show=FALSE,muinit=c(5,2),siginit=c(2,-0.1),shinit=c(0,2))
  
  checkEqualsNumeric(ismod$mle,coef(mod),tolerance = tol,msg="gev: covariates in mu phi and xi: point ests")
  checkEqualsNumeric(ismod$se,sqrt(diag(mod$cov)),tolerance = tol,msg="gev: covariates in mu phi and xi: std errs")
  
  ####################################################################
  # 5.5 GEV:  Check that using priors gives expected behaviour when covariates are included.
  
  # Tests for xi being drawn to 0
  
  myb <- rep(c(-0.1,0.1),each=5)
  data <- makeDataGev(u=0,a=1,b=-0.5+myb,n=3000)
  data$b <- myb
  
  gp1 <- list(c(0, 0, 0, 0), diag(c(10^4, 10^4, 0.25, 0.25)))
  gp2 <- list(c(0, 0, 0, 0), diag(c(10^4, 10^4, 0.25, 0.01)))
  
  mod0 <- evm(y,family=gev,data=data,xi=~b,penalty="none")
  mod1 <- evm(y,family=gev,data=data,xi=~b,priorParameters=gp1)
  mod2 <- evm(y,family=gev,data=data,xi=~b,priorParameters=gp2)
  
  checkTrue(all(abs(coef(mod0)[3:4]) > abs(coef(mod1)[3:4])),
            msg="gev: with covariates, xi drawn to zero")
  checkTrue(abs(coef(mod1)[4]) > abs(coef(mod2)[4]),
            msg="gev: with covariates, xi drawn to zero")
  
  # Tests for phi being drawn to 0
  
  mya <- seq(0.1,1,len=10)
  data <- makeDataGev(u=0,a=-3 + mya,b=-0.1,n=3000)
  data$a <- rep(mya, len=nrow(data))
  
  gp4 <- list(c(0, 0, 0, 0), diag(c(10^4, 1, 1, 10^4)))
  gp5 <- list(c(0, 0, 0, 0), diag(c(10^4, 0.1, 0.1, 10^4)))
  
  mod3 <- evm(y,family=gev,data=data,phi=~a,penalty="none")
  mod4 <- evm(y,family=gev,data=data,phi=~a,priorParameters=gp4)
  mod5 <- evm(y,family=gev,data=data,phi=~a,priorParameters=gp5)
  
  checkTrue(all(abs(coef(mod3)[2:3]) > abs(coef(mod4)[2:3])),
            msg="gev: with covariates, phi being drawn to 0")
  checkTrue(all(abs(coef(mod4)[2:3]) > abs(coef(mod5)[2:3])),
            msg="gev: with covariates, phi being drawn to 0")
  
  # Tests for xi being drawn to 2
  myb <- rep(c(-0.1,0.1),each=5)
  data <- makeDataGev(u=0,a=1,b=myb,n=3000)
  
  gp7 <- list(c(0, 0, 1, 1), diag(c(10^4, 10^4, 0.25, 0.25)))
  gp8 <- list(c(0, 0, 1, 1), diag(c(10^4, 10^4, 0.05, 0.05)))
  
  mod6 <- evm(y,family=gev,data=data,xi=~b,penalty="none")
  mod7 <- evm(y,family=gev,data=data,xi=~b,priorParameters=gp7)
  mod8 <- evm(y,family=gev,data=data,xi=~b,priorParameters=gp8)
  
  checkTrue(all(abs(1 - coef(mod6)[3:4]) > abs(1 - coef(mod7)[3:4])),
            msg="gev: with covariates, xi drawn to 1")
  checkTrue(all(abs(1 - coef(mod7)[3:4]) > abs(1 - coef(mod8)[3:4])),
            msg="gev: with covariates, xi drawn to 1")
  
  # Tests for mu being drawn to 4
  
  myu <- seq(0.1,1,len=10)
  data <- makeDataGev(2 + myu,a=1,b=-0.1,n=3000)
  data$u <- rep(myu, len=nrow(data))
  
  gp10 <- list(c(4, 0, 0, 0), diag(c(1,   10^4, 10^4, 10^4)))
  gp11 <- list(c(4, 0, 0, 0), diag(c(0.1, 10^4, 10^4, 10^4)))
  
  mod9 <-  evm(y,family=gev,data=data,mu=~u,penalty="none")
  mod10 <- evm(y,family=gev,data=data,mu=~u,priorParameters=gp10)
  mod11 <- evm(y,family=gev,data=data,mu=~u,priorParameters=gp11)
  
  checkTrue(abs(4 - coef(mod9)[1])  > abs(4 - coef(mod10)[1]),
            msg="gev: with covariates, mu drawn to 4")
  checkTrue(abs(4 - coef(mod10)[1]) > abs(4 - coef(mod11)[1]),
            msg="gev: with covariates, mu drawn to 4")
  
  
}

test.evmBoot <- function(){
  tol <- 0.1
  
  for(Family in list(gpd,gev)){
    set.seed(20130615)
    
    pst <- function(msg) texmexPst(msg,Family=Family)
    
    u    <- switch(Family$name,GPD=30,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    fit <- evm(data,th=u,family=Family, penalty="none")
    boot <- evmBoot(fit, R=200, trace=1000)
    co <- coef(fit)
    rep <- boot$replicates
    scaleColumn <- switch(Family$name,GPD=1,GEV=2)
    rep[,scaleColumn] <- exp(rep[, scaleColumn])
    
    # Compare bootstrap standard errors with those given by Coles
    # pages 59 and 85 for GEV and GPD resp
    
    bse <- apply(rep, 2, sd)
    cse <- switch(Family$name,GPD=c(.958432, .101151),GEV=c(0.02792848, 0.02024846, 0.09823441))
    
    checkEqualsNumeric(cse,bse,tolerance=tol,
                       msg=pst("evmBoot: bootstrap se of parameter estimates matches Coles"))
    
    ## Check penalization works - set harsh penalty and do similar
    ## checks to above
    
    pp <- switch(Family$name,GPD=list(c(0, .5), diag(c(.5, .05))), GEV=list(c(5, 0, .5), diag(c(.5, .5, .05))))
    fit <- evm(data, th=u, penalty="none", priorParameters=pp,family=Family)
    boot <- evmBoot(fit, R=1000, trace=1100)
    
    bse <- apply(boot$replicates, 2, sd)
    rse <- bse / fit$se
    rse <- ifelse(rse < 1, 1/rse, rse)
    checkTrue(max(rse) < 1.1, msg=pst("evmBoot: SEs with xi in model, with penalty applied"))
    
    best <- apply(boot$replicates, 2, median)
    fest <- coef(fit)
    rdiff <- abs((best - fest)/fest)
    checkTrue(all(rdiff < 0.06), msg=pst("evmBoot: medians in line with point ests, with penalty applied"))
    
    ##################################################################
    # models with covariates. Due to apparent instability
    # of the Hessian in some cases, allow some leeway
    
    n <- 1000
    mu <- 1
    phi <- 5
    xi <- 0.05
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    
    test <- function(boot,fit,txt){
      bse <- apply(boot$replicates, 2, sd)
      rse <- bse / fit$se
      rse <- ifelse(rse < 1, 1/rse, rse)
      checkTrue(max(rse) < 1.5, msg=pst(paste("evmBoot: SEs with covariates in",txt)))
      
      best <- apply(boot$replicates, 2, median)
      fest <- coef(fit)
      rdiff <- abs((best - fest)/fest)
      
      checkTrue(all(rdiff < 0.2), msg=pst(paste("evmBoot: medians in line with point ests, covariates in",txt)))
    }
    
    param <- switch(Family$name,GPD=cbind(2+X[,1],xi),GEV=cbind(mu,2+X[,1],xi))
    start <- switch(Family$name,GPD=c(2,1,xi),GEV=c(mu,2,1,xi))
    X$Y <- Family$rng(n,param,list(threshold=th))
    
    fit <- evm(Y,data=X,phi=~a,th=th,family=Family,start=start)
    boot <- evmBoot(fit, R=200, trace=201)
    test(boot,fit,"phi")
    
    param <- switch(Family$name,GPD=cbind(phi,0.1+X[,2]),GEV=cbind(mu,phi,0.1+X[,2]))
    start <- switch(Family$name,GPD=c(phi,0.1,1),GEV=c(mu,phi,0.1,1))
    X$Y <- Family$rng(n,param,list(threshold=th))
    
    fit <- evm(Y,data=X,xi=~b,th=th,family=Family,start=start)
    boot <- evmBoot(fit, R=200, trace=201)
    test(boot,fit,"xi")
  }
}

test.extremalIndex <- function(){
  tol <- 0.0001
  th <- quantile(rain,seq(0.7,0.99,len=10))
  for(i in 1:length(th)){
    texmex.ei <- extremalIndex(rain,threshold=th[i])
    Ferro.ei  <- texmex:::.extRemes.exi.intervals(rain > th[i])
    
    Ferro.clust <- texmex:::.extRemes.decluster.intervals(rain> th[i], Ferro.ei)
    texmex.clust <- declust(texmex.ei)
    
    Ferro.runs <-  texmex:::.extRemes.decluster.runs(rain> th[i], 3)
    texmex.runs <- declust(rain,threshold=th[i],r=3,verbose=FALSE)
    
    checkEqualsNumeric(texmex.ei$EIintervals, Ferro.ei,
                       tolerance = tol,msg="extremalIndex: extRemes implementation")
    checkEqualsNumeric(texmex.clust$sizes, Ferro.clust$size,
                       tolerance = tol,msg="extremalIndex: declustering")
    
    checkEqualsNumeric(texmex.runs$nCluster,Ferro.runs$nc,msg="extremalIndex: runs declustering nc")
    checkEqualsNumeric(texmex.runs$sizes,Ferro.runs$size,msg="extremalIndex: runs declustering sizes")
  }
  
  # check passing data through data frames
  
  data <- data.frame(RAIN=rain[1:1000], rnorm=rnorm(1000), num=1:1000)
  extremalIndexRangeFit(RAIN, data,verbose=FALSE,nboot=20,nint=7)
  extremalIndexRangeFit(data$RAIN,verbose=FALSE,nboot=20,nint=7)
  
  data.de <- declust(RAIN,data=data,th=th[1],verb=FALSE)
  resp.de <- declust(data$RAIN,th=th[1],verb=FALSE)
  
  data.ei <- extremalIndex(RAIN,data=data,threshold=th[1])
  resp.ei <- extremalIndex(data$RAIN,threshold=th[1])
  
  checkEqualsNumeric(data.ei$EIintervals,resp.ei$EIintervals,tolerance=tol,msg="extremalIndex: using data frame to pass response")
  checkEqualsNumeric(data.de$clusters,resp.de$clusters,tolerance=tol,msg="extremalIndex: using data frame to pass numeric response to declustering")
  
  # test covariate fitting
  
  ei <- extremalIndex(SO2,data=winter,threshold=20)
  d <- declust(ei)
  evm(d,phi=~NO)
  
  checkEqualsNumeric(662.9508, AIC(evm(d,phi=~NO)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(662.8874, AIC(evm(d,phi=~NO2)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(651.8747, AIC(evm(d,phi=~O3)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(663.0015, AIC(evm(d,phi=~PM10)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(651.7874, AIC(evm(d,phi=~O3,xi=~NO)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(653.2512, AIC(evm(d,phi=~O3,xi=~NO2)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(653.6385, AIC(evm(d,phi=~O3,xi=~O3)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(652.9238, AIC(evm(d,phi=~O3,xi=~PM10)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  
}

test.gpd.info <- function(){
  lmod <- evm(ALT.M, data=liver, qu=.5, xi=~I(240*as.numeric(dose)),
              cov="numeric")
  checkTrue(all(sqrt(diag(solve(gpd.info(lmod)))) > 0),
            msg="gpd.inf: SDs positive")
  
  # Check equality to numerical approximation in big samples
  set.seed(20110923)
  tol <- 10^(-3)
  for (i in 1:10){
    x <- rt(10000, 10)
    junk <- evm(x, qu=.9, penalty="none", cov="numeric")
    msg <- paste("gpd.info: t", i, "equality to numerical", sep="")
    checkEqualsNumeric(junk$cov, solve(gpd.info(junk)), tolerance=tol,
                       msg=msg)
    
    # check estimation when we have a penalty
    gp1 <- list(c(0, 0), diag(c(10^4, .05)))
    gp2 <- list(c(0, 0), diag(c(.1, 10^4)))
    junk1 <- evm(x, qu=.9, priorParameters = gp1, cov="numeric")
    junk2 <- evm(x, qu=.9, priorParameters = gp2, cov="numeric")
    msg1 <- paste("gpd.info: t", i, "equality to numerical, penalty on xi", sep="")
    msg2 <- paste("gpd.info: t", i, "equality to numerical, penalty on phi", sep="")
    tol <- 0.01
    checkEqualsNumeric(junk1$cov, solve(gpd.info(junk1)), tolerance=tol, msg=msg1)
    checkEqualsNumeric(junk2$cov, solve(gpd.info(junk2)),
                       tolerance=tol, msg=msg2)
    
    # check estimation when we have covariates
    n <- 10000
    x <- 1/runif(n)
    data <- data.frame(x=x,y=rexp(n,exp(2 + x)))
    
    junk3 <- evm(y,data=data,phi =~ x,th=0)
    msg3 <- paste("gpd.info: t",i,"equality to numerical, covariates in phi",sep="")
    checkEqualsNumeric(junk3$cov, solve(gpd.info(junk3)), tolerance=tol, msg=msg3)
    
    x <- runif(n,-0.5,0.5)
    data <- data.frame(x=x,y = rgpd(n,sigma = exp(3+2*x), xi=x))
    
    junk4 <- evm(y,data=data,phi=~x, xi = ~ x,th=0)
    msg4 <- paste("gpd.info: t",i,"equality to numerical, covariates in phi and xi",
                  sep="")
    checkEqualsNumeric(junk4$cov, solve(gpd.info(junk4)), tolerance=tol, msg=msg4)
  }
}

test.gpdRangeFit <- function(){
  par(mfrow=c(2,1))
  res <- gpdRangeFit(rain, umin=0, umax=50, nint=20)
  plot(res, pch=16, main=c("Figure 4.2 of Coles (2001)",""), addNexcesses=FALSE)
  checkEquals(names(res),c("th", "par","hi","lo","data"),msg=pst("gpdRangeFit: returned list with correct names"))
}

test.mexDependence <- function(){
  smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
  wmarmod <- migpd(winter, mqu=.7,  penalty="none")
  
  mySdepO3 <- mexDependence(smarmod,which=1,dqu=0.7,margins="gumbel",constrain=FALSE)
  myWdepO3 <- mexDependence(wmarmod,which=1,dqu=0.7,margins="gumbel",constrain=FALSE)
  
  mySdepNO2 <- mexDependence(smarmod,which=2,dqu=0.7,margins="gumbel",constrain=FALSE)
  myWdepNO2 <- mexDependence(wmarmod,which=2,dqu=0.7,margins="gumbel",constrain=FALSE)
  
  mySdepNO <- mexDependence(smarmod,which=3,dqu=0.7,margins="gumbel",constrain=FALSE)
  myWdepNO <- mexDependence(wmarmod,which=3,dqu=0.7,margins="gumbel",constrain=FALSE)
  
  mySdepSO2 <- mexDependence(smarmod,which=4,dqu=0.7,margins="gumbel",constrain=FALSE)
  myWdepSO2 <- mexDependence(wmarmod,which=4,dqu=0.7,margins="gumbel",constrain=FALSE)
  
  mySdepPM10 <- mexDependence(smarmod,which=5,dqu=0.7,margins="gumbel",constrain=FALSE)
  myWdepPM10 <- mexDependence(wmarmod,which=5,dqu=0.7,margins="gumbel",constrain=FALSE)
  
  
  jhSdepO3 <- matrix(c(
    0.56627103,  0.37272912, 0.0000000, 0.0000000,
    0.22029334,  0.36865296, 0.0000000, 0.0000000,
    0.28193999, -0.26731823, 0.0000000, 0.0000000,
    0.46293139, -0.23387868, 0.0000000, 0.0000000),byrow=FALSE,nrow=4)
  
  jhSdepNO2 <- matrix(c(
    0.49290567,  0.22236302, 0.0000000, 0.0000000,
    0.38571246,  0.34379705, 0.0000000, 0.0000000,
    0.22026515, -0.17494068, 0.0000000, 0.0000000,
    0.45455612,  0.22411795, 0.0000000, 0.0000000),byrow=FALSE,nrow=4)
  
  jhSdepNO <- matrix(c(
    0.43149222,  0.34033851, 0.0000000, 0.0000000,
    0.49992799,  0.21878814, 0.0000000, 0.0000000,
    0.19724402,  0.23839660, 0.0000000, 0.0000000,
    0.50384850,  0.18227312, 0.0000000, 0.0000000),byrow=FALSE,nrow=4)
  
  jhSdepSO2 <- matrix(c(
    0.24400046, -0.02162792, 0.0000000, 0.0000000,
    0.08769596, -0.14758165, 0.0000000, 0.0000000,
    0.00000000, -0.04461209, 0.6865857, 0.4201682,
    0.35364948,  0.02338747, 0.0000000, 0.0000000), byrow=FALSE,nrow=4)
  
  jhSdepPM10 <- matrix(c(
    0.08302144,  0.16604598, 0.0000000, 0.0000000,
    0.00000000,  0.57387887, 0.0000000, 0.0000000,
    0.15208086,  0.32264497, 0.0000000, 0.0000000,
    0.00000000,  0.43255493, 0.0000000, 0.0000000), byrow=FALSE,nrow=4)
  
  
  
  jhWdepO3 <- matrix(c(
    0.00000000,  0.008072046,  0.00000000, 0.0000000,
    0.00000000,  0.034283871,  0.00000000, 0.0000000,
    0.00000000, -0.188517544,  5.14775893, 1.0000000,
    0.00000000, -0.026874734,  0.05011460, 0.1075632),byrow=FALSE,nrow=4)
  
  jhWdepNO2 <- matrix(c(
    0.00000000, -0.553608371, -0.06047238, 0.4967213,
    0.81920276,  0.529272235,  0.00000000, 0.0000000,
    0.32246150,  0.335335739,  0.00000000, 0.0000000,
    0.85746906,  0.085265792,  0.00000000, 0.0000000),byrow=FALSE,nrow=4)
  
  jhWdepNO <- matrix(c(
    0.00000000, -0.504344703, -1.41890419, 0.0000000,
    0.75819233,  0.378119827,  0.00000000, 0.0000000,
    0.32199902, -0.350339706,  0.00000000, 0.0000000,
    0.73227271, -0.105822435,  0.00000000, 0.0000000),byrow=FALSE,nrow=4)
  
  jhWdepSO2 <- matrix(c(
    0.00000000, -0.485253436, -1.27253412, 0.0000000,
    0.00000000, -0.018577905,  0.63501876, 0.3862878,
    0.00000000,  0.000000000,  0.76856266, 0.4916768,
    0.03626605, -0.316472032,  0.00000000, 0.0000000),byrow=FALSE,nrow=4)
  
  jhWdepPM10 <- matrix(c(
    0.00000000,  0.064075145,  0.00000000, 0.0000000,
    0.86288696,  0.584629421,  0.00000000, 0.0000000,
    0.59510081,  0.569002154,  0.00000000, 0.0000000,
    0.10412199,  0.207529741,  0.00000000, 0.0000000),byrow=FALSE,nrow=4)
  
  tol <- 0.23
  if(FALSE){
    par(mfrow=c(2,5))
    plot(jhWdepO3,  myWdepO3$dependence$coefficients);abline(0,1)
    plot(jhWdepNO2, myWdepNO2$dependence$coefficients);abline(0,1)
    plot(jhWdepNO,  myWdepNO$dependence$coefficients);abline(0,1)
    plot(jhWdepSO2, myWdepSO2$dependence$coefficients);abline(0,1)
    plot(jhWdepPM10,myWdepPM10$dependence$coefficients);abline(0,1)
    
    plot(jhSdepO3,  mySdepO3$dependence$coefficients);abline(0,1)
    plot(jhSdepNO2, mySdepNO2$dependence$coefficients);abline(0,1)
    plot(jhSdepNO,  mySdepNO$dependence$coefficients);abline(0,1)
    plot(jhSdepSO2, mySdepSO2$dependence$coefficients);abline(0,1)
    plot(jhSdepPM10,mySdepPM10$dependence$coefficients);abline(0,1)
  }
  
  checkEqualsNumeric(jhWdepO3,  myWdepO3$dependence$coefficients[1:4,],  tolerance=tol,msg="mexDependence: Winter O3")
  checkEqualsNumeric(jhWdepNO2, myWdepNO2$dependence$coefficients[1:4,], tolerance=tol,msg="mexDependence: Winter NO2")
  checkEqualsNumeric(jhWdepNO,  myWdepNO$dependence$coefficients[1:4,],  tolerance=tol,msg="mexDependence: Winter NO")
  checkEqualsNumeric(jhWdepSO2, myWdepSO2$dependence$coefficients[1:4,], tolerance=tol,msg="mexDependence: Winter SO2")
  checkEqualsNumeric(jhWdepPM10,myWdepPM10$dependence$coefficients[1:4,],tolerance=tol,msg="mexDependence: Winter PM10")
  
  checkEqualsNumeric(jhSdepO3,  mySdepO3$dependence$coefficients[1:4,],  tolerance=tol,msg="mexDependence: Summer O3")
  checkEqualsNumeric(jhSdepNO2, mySdepNO2$dependence$coefficients[1:4,], tolerance=tol,msg="mexDependence: Summer NO2")
  checkEqualsNumeric(jhSdepNO,  mySdepNO$dependence$coefficients[1:4,],  tolerance=tol,msg="mexDependence: Summer NO")
  checkEqualsNumeric(jhSdepSO2, mySdepSO2$dependence$coefficients[1:4,], tolerance=tol,msg="mexDependence: Summer SO2")
  checkEqualsNumeric(jhSdepPM10,mySdepPM10$dependence$coefficients[1:4,], tolerance=tol,msg="mexDependence: Summer PM10")
  
  # test functionality with 2-d data
  
  wavesurge.fit <- migpd(wavesurge,mqu=.8)
  dqu <- 0.8
  which <- 1
  wavesurge.mex <- mexDependence(wavesurge.fit,which=which,dqu=dqu)
  op <- options(warn=-1)
  wavesurge.dep <- mexDependence(wavesurge.fit,which=which)
  options(op)
  
  checkEquals(wavesurge.mex[1:2],wavesurge.dep[1:2],msg="mexDependence: missing dqu argument")
  checkEqualsNumeric(dim(wavesurge.mex$dependence$Z),c(578,1),msg="mexDependence: execution for 2-d data")
  checkEqualsNumeric(wavesurge.mex$dependence$dqu, dqu, msg="mexDependence: execution for 2-d data")
  checkEqualsNumeric(wavesurge.mex$dependence$which,which,msg="mexDependence: execution for 2-d data")
  
  # test specification of starting values
  
  which <- 2
  dqu <- 0.8
  summer.fit <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
  summer.mex1 <- mexDependence(summer.fit,which=which,dqu=dqu)
  summer.mex2 <- mexDependence(summer.fit,which=which,dqu=dqu,start=summer.mex1)
  summer.mex3 <- mexDependence(summer.fit,which=which,dqu=dqu,start=summer.mex1$dependence$coefficients[1:2,])
  summer.mex4 <- mexDependence(summer.fit,which=which,dqu=dqu,start=c(0.01,0.03))
  
  tol <- 0.005
  checkEqualsNumeric(summer.mex1$dependence$coefficients,summer.mex2$dependence$coefficients,tolerance=tol,msg="mexDependence: summer starting values: class(start)='mex'")
  checkEqualsNumeric(summer.mex1$dependence$coefficients,summer.mex3$dependence$coefficients,tolerance=tol,msg="mexDependence: summer starting values: start= matrix ")
  checkEqualsNumeric(summer.mex1$dependence$coefficients,summer.mex4$dependence$coefficients,tolerance=tol,msg="mexDependence: summer starting values: start= vector")
  
  # test laplace constrained estimation against independent coding of implementation by Yiannis Papastathopoulos (see test.functions.R file)
  # do all the possible pairs generated by the winter data set.  This covers negative and positive dependence with cases on and off the boundary.
  
  set.seed(20111010)
  pairs <- expand.grid(1:5,1:5)
  pairs <- as.matrix(pairs[pairs[,1] != pairs[,2], 2:1])
  margins <- casefold("laplace")
  marTransform <- c("mixture","empirical")
  Dthresh <- c(1,2,2.5)
  Mquant <- 0.7
  v <- 10
  k <- 5
  
  for(marTrans in marTransform){
    for(dth in Dthresh){
      HTSest <- KTPest <- matrix(0,ncol=2,nrow=dim(pairs)[1])
      
      for(i in 1:dim(pairs)[1]){
        mar.model <- migpd(winter[,pairs[i,]],mqu=Mquant,penalty="none")
        mar.trans <- mexTransform(mar.model,margins=margins,method=marTrans)
        X <- list(mar.trans$transformed)
        Dqu <- 1 - mean(mar.trans$transformed[,1] > dth)
        init <- initial_posneg(D=1,listdata=X,u=dth,v=v)
        a <- estimate_HT_KPT_joint_posneg_nm(pars=init,x=dth,listr=X,params=FALSE,k=k,v=v)
        KTPest[i,] <- a$par
        b <- mexDependence(mar.trans,which=1,dqu=Dqu,margins=margins,constrain=TRUE,v=v,
                           marTransform=marTrans,start=init,nOptim=k)
        HTSest[i,] <- coef(b)$dependence[1:2]
      }
      
      checkEqualsNumeric(KTPest,HTSest,tolerance=tol,mesg=paste("mexDependence: constrained estimation threshold",i))
    }
  }
}

test.mexRangeFit <- function(){
  
  which <- 2
  quantiles <- seq(0.5,0.9,length=5)
  wmarmod <- migpd(winter, mqu=.7,  penalty="none")
  wmexmod.gum <- mexDependence(wmarmod, dqu=quantiles[1], margins="gumbel", constrain=FALSE,which=which)
  wmexmod.lap <- mexDependence(wmarmod, dqu=quantiles[1], margins="laplace",v=5,which=which)
  
  R <- 3
  mrf1 <- mexRangeFit(wmarmod,quantiles = quantiles,which=which,R=R,trace=R+1,v=5)
  mrf2 <- mexRangeFit(wmexmod.gum,quantiles = quantiles,R=R,trace=R+1)
  mrf3 <- mexRangeFit(wmexmod.lap,quantiles = quantiles,R=R,trace=R+1)
  
  checkException(mexRangeFit(TRUE,which=2),silent=TRUE,msg="mexRangeFit: exception handle")
  checkException(mexRangeFit(5,which=1),silent=TRUE,msg="mexRangeFit: exception handle")
  
  checkEquals(mrf1$ests[[1]][1:2],wmexmod.lap[1:2])
  checkEquals(mrf2$ests[[1]][1:2],wmexmod.gum[1:2])
  checkEquals(mrf3$ests[[1]][1:2],wmexmod.lap[1:2])
  
  # now 2-d data
  
  mqu <- .7
  wavesurge.fit <- migpd(wavesurge,mqu=mqu)
  m <- mexDependence(wavesurge.fit,which=1,dqu=mqu)
  mrf4 <- mexRangeFit(wavesurge.fit,which=1,margins="laplace",R=R,trace=R+1)
  mrf5 <- mexRangeFit(m,R=R,trace=R+1)
  checkEquals(mrf4$ests[[2]][1:2],mrf5$ests[[2]][1:2])
  
  # test specification of starting values
  R <- 5
  qu <- c(0.5,0.7,0.9)
  mrf6 <- mexRangeFit(wavesurge.fit,which=1,margins="laplace",constrain=TRUE, start=c(0.01,0.01),R=R,trace=R+1,quantiles = qu)
  mrf7 <- mexRangeFit(wavesurge.fit,which=2,margins="laplace",constrain=TRUE, start=m,R=R,trace=R+1,quantiles = qu)
  
  # test plotting
  
  par(mfrow=c(2,2))
  plot(mrf6,main="start=(0.01,0.01)",addNexcesses=FALSE)
  plot(mrf7,main=paste("start=",signif(coef(m)$dependence[1:2],2)),addNexcesses=FALSE)
}

test.migpd <- function(){
  
  # values from Heffernan and Tawn (2004) Table 4.
  # Note values in published Table 4 for u_{Xi} in cols NO2 and NO Winter were reversed.
  
  htsummer <- rbind(mqu=c(43, 43, 66.1, 22, 46),
                    mth = c(.9, .7, .7, .85, .7),
                    sigma = c(15.8, 9.1, 32.2, 42.9, 22.8),
                    xi = c(-.29, .01, .02, .08, .02))
  
  htwinter <- rbind(mqu=c(28, 49, 151.6, 23, 53),
                    mth = rep(.7, 5),
                    sigma = c(6.2, 9.3, 117.4, 19.7, 37.5),
                    xi = c(-.37, -.03, -.09, .11, -.2))
  
  summer.gpd <- summary(migpd(summer, mqu=htsummer[2,],penalty="none"),verbose=FALSE)
  winter.gpd <- summary(migpd(winter, mqu=htwinter[2,],penalty="none"),verbose=FALSE)
  
  tol <- c(1, 0.05, .5, 0.5)
  for(i in 1:4){
    checkEqualsNumeric(htsummer[i,], summer.gpd[i,],tolerance=tol[i],msg=paste("migpd: Table 4 summer",i))
    checkEqualsNumeric(htwinter[i,], winter.gpd[i,],tolerance=tol[i],msg=paste("migpd: Table 4 winter",i))
  }
  
  # check excecution for 2-d data
  
  wavesurge.fit <- migpd(wavesurge,mqu=.7)
  checkEqualsNumeric(wavesurge.fit$models$wave$loglik, evm(wavesurge$wave,qu=0.7)$loglik,
                     tolerance=0.001,msg="migpd: 2-d data gpd fit wave")
  
}

test.migpdCoefs <- function(){
  
  liver <- liver
  liver$ndose <- as.numeric(liver$dose)
  require(MASS,quietly=TRUE) # For rlm
  
  ralt <- resid(rlm(log(ALT.M) ~ log(ALT.B) + ndose, data=liver))
  rast <- resid(rlm(log(AST.M) ~ log(AST.B) + ndose, data=liver))
  ralp <- resid(rlm(log(ALP.M) ~ log(ALP.B) + ndose, data=liver))
  rtbl <- resid(rlm(log(TBL.M) ~ log(TBL.B) + ndose, data=liver))
  
  rliver <- data.frame(alt=ralt, ast=rast, alp=ralp, tbl=rtbl, ndose=liver$ndose)
  
  Dmod <- migpd(rliver[rliver$ndose == 4, 1:4], mqu=.7) # Model for dose D
  
  oldALTco <- coef(Dmod)[3:4, 1]
  
  altgpd <- evm(alt, qu=.7, xi = ~ ndose, data=rliver)
  astgpd <- evm(ast, qu=.7, xi = ~ ndose, data=rliver)
  
  altco <- c(coef(altgpd)[1], coef(altgpd)[2] + 4 * coef(altgpd)[3])
  astco <- c(coef(astgpd)[1], coef(astgpd)[2] + 4 * coef(astgpd)[3])
  
  # Change one set of coefficients
  lmod <- migpdCoefs(Dmod, which="alt", list(altco))
  
  newALTco <- coef(lmod)[3:4, 1]
  newALTco[1] <- log(newALTco[1]) # For comparison with altco
  oldALTco[1] <- log(oldALTco[1])
  
  checkEqualsNumeric(altco, newALTco, msg="migpdCoefs: change one set of coefficients")
  checkTrue(all(newALTco != oldALTco), msg="migpdCoefs: change one set of coefficients")
  
  # Change 2 sets of coefficients at once
  
  lmod <- migpdCoefs(Dmod, which=c("alt", "ast"), coefs=list(altco, astco))
  
  newCo <- coef(lmod)[3:4, 1:2]
  oldCo <- coef(Dmod)[3:4, 1:2]
  
  newCo[1,] <- log(newCo[1,])
  oldCo[1,] <- log(oldCo[1,])
  
  checkEqualsNumeric(cbind(altco, astco),newCo, msg="migpdCoefs: change two set of coefficients at once")
  checkTrue(all(newCo != oldCo), "migpdCoefs: change two set of coefficients at once")
}

test.mrl <- function(){
  par(mfrow=c(1,1))
  res <- mrl(rain)
  res <- plot(res, , main="Figure 4.1 of Coles (2001)")
  checkEquals(res, NULL, msg="mrlPlot: check execution")
}

test.MCS <- function(){
  myMCS <- function(x,p){
    # First and second args are
    # x (dxn matrix) and p (vector of probabilities).
    
    n <- dim(x)[2]
    d <- dim(x)[1]
    u <- t(apply(x,1,rank))/(n+1)# schmid and schmidt use n not n+1
    
    rho <- numeric(length(p))
    for(k in 1:length(p)) {
      Diff <- p[k] - u
      Diff[Diff<0] <- 0
      Prod <- apply(Diff,2,prod)
      
      num <- mean(Prod) - ((p[k]^2)/2)^d
      den <- (p[k]^(d+1)) / (d+1) - (0.5*p[k]^2)^d
      rho[k] <- num/den
    }
    
    rho
  }
  # simulated data - dimension 2
  n <- 1000
  by <- 0.01
  p <- seq(by,1-by,by=by)
  data <- rbind(rnorm(n),rnorm(n))
  
  tmRl <- MCS(t(data),p)
  myRl <- myMCS(data,tmRl$p)
  
  checkEqualsNumeric(myRl,tmRl$mcs,msg="MCS: independent normal data")
  checkEqualsNumeric(p,tmRl$p, msg="MCS: mathching p argument")
  
  # winter air pollution data - dimension 5
  tmWinterMCS <- MCS(winter,p)
  myWinterMCS <- myMCS(t(winter),p)
  checkEqualsNumeric(myWinterMCS,tmWinterMCS$mcs,msg="MCS: winter air pollution data")
  
  # summer airpollution data - dimension 5
  tmSummerMCS <- MCS(summer,p)
  mySummerMCS <- myMCS(t(summer),p)
  checkEqualsNumeric(mySummerMCS,tmSummerMCS$mcs,msg="MCS: summer air pollution data")
}

test.pgev <- function() {
  probabilities <- runif(10)
  
  xi.values <- c(0, seq(-5, 5, length.out=10))
  
  core.sanity.test <- function(xi) {
    randoms <- sort(rgev(10, 0, 1, xi))
    checkTrue(all(diff(pgev(randoms, 0, 1, xi)) >= 0),
              "pgev: ascending")
    checkTrue(all(diff(pgev(randoms, 0, 1, xi, lower.tail=FALSE)) <= 0),
              "pgev: descending")
    
    checkEqualsNumeric(log(pgev(randoms, 0, 1, xi)),
                       pgev(randoms, 0, 1, xi, log.p=TRUE),
                       "pgev: log.p (1)")
    checkEqualsNumeric(log(pgev(randoms, 0, 1, xi, lower.tail=FALSE)),
                       pgev(randoms, 0, 1, xi, lower.tail=FALSE,
                            log.p=TRUE),
                       "pgev: log.p (2)")
    
    mu <- runif(1, -5, 5)
    sigma <- rexp(1)
    
    checkEqualsNumeric(pgev(randoms, 0, 1, xi),
                       pgev(mu + sigma * randoms, mu, sigma, xi),
                       "pgev: shift and scale")
  }
  
  qgev.comparison <- function(xi) {
    quantiles <- qgev(probabilities, 0, 1, xi)
    my.probs  <- pgev(quantiles, 0, 1, xi)
    checkEqualsNumeric(my.probs, probabilities, "pgev: straight test")
    
    quantiles <- qgev(probabilities, 0, 1, xi, lower.tail=FALSE)
    my.probs  <- pgev(quantiles, 0, 1, xi, lower.tail=FALSE)
    checkEqualsNumeric(my.probs, probabilities, "pgev: lower tail")
    
    my.probs.2 <- pgev(quantiles, 0, 1, xi, lower.tail=TRUE)
    checkEqualsNumeric(probabilities + my.probs.2,
                       rep(1, length(probabilities)),
                       "pgev: tail flip")
  }
  
  lapply(xi.values, core.sanity.test)
  lapply(xi.values, qgev.comparison)
}
test.pgpd <- function(){
  
  evd.pgpd <- texmex:::.evd.pgpd
  myTest <- function(sig,xi,thresh,msg){
    myp <- sapply(1:nreps,function(i) pgpd(x[,i], sig[i], xi[i],u=thresh[i]))
    ep <- sapply(1:nreps, function(i) evd.pgpd(x[,i], loc=thresh[i], scale=sig[i], shape=xi[i]))
    checkEqualsNumeric(ep,myp,msg=msg)
  }
  
  set.seed(20101111)
  
  #*************************************************************
  # 6.7. Test pgpd. Note that .evd.pgpd is NOT vectorized.
  
  nreps <- 100
  nsim <- 1000
  p <- matrix(runif(2*nreps, -1, 1),ncol=2)
  p[, 1] <- p[, 1] + 1
  thresh <- rep(0,nreps)
  
  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))
  
  myTest(sig=p[,1], xi=p[,2],thresh=thresh, msg="pgpd: random xi")
  
  #*************************************************************
  # 6.8. Test pgpd when some or all of xi == 0
  
  p[sample(1:nreps,nreps/2),2] <- 0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="pgpd: some zero xi")
  
  p[,2] <-  0
  x <- sapply(1:nreps,function(i)rgpd(nsim,sigma=p[i,1],xi=p[i,2],u=thresh[i]))
  myTest(sig=p[,1], xi=p[,2], thresh=thresh, msg="pgpd: all zero xi")
  
  #*************************************************************
  # 6.9. Test vectorization of pgpd.
  
  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)
  thresh <- rnorm(nsim)
  
  x <- rgpd(nsim, sig, xi,u=thresh)
  myp <- pgpd(x, sig, xi,u=thresh)
  
  ep <- sapply(1:nsim, function(i)evd.pgpd(x[i], loc=thresh[i], scale=sig[i], shape=xi[i]))
  checkEqualsNumeric(ep,myp,msg="pgpd: vectorisation")
  
  #*************************************************************
  # 6.10 test log.p argument
  
  lp <- pgpd(x,sig,xi,u=thresh,log.p=TRUE)
  checkEqualsNumeric(myp,exp(lp),msg="pgpd: log probabilities")
  
  #*************************************************************
  # 6.11 test lower tail argument
  
  sp <- pgpd(x,sig,xi,u=thresh,lower.tail=FALSE)
  checkEqualsNumeric(myp,1-sp,msg="pgpd: lower tail")
  
  ## check pgpd when q < threshold
  upperProb <- pgpd(0, 1, 1, u=0.5, lower.tail=TRUE)
  checkEqualsNumeric(upperProb, 0, msg="pgpd: value below threshold (1)")
  
  lowerProb <- pgpd(0, 1, 1, u=0.5, lower.tail=FALSE)
  checkEqualsNumeric(upperProb, 0, msg="pgpd: value below threshold (2)")
  
  ## check pgpd when xi < 0 and value above upper limit
  
  xi <- -2.3
  upperProb <- pgpd(-2/xi, 1, xi, u=0, lower.tail=TRUE)
  checkEqualsNumeric(upperProb, 1, msg="pgpd: negative xi (1)")
  
  lowerProb <- pgpd(-2/xi, 1, xi, u=0, lower.tail=FALSE)
  checkEqualsNumeric(lowerProb, 0, msg="pgpd: negative xi (2)")
}

test.plot.bootmex <- function(){
  
  set.seed(3141593)
  
  # 2-d wavesurge data
  
  wavesurge.fit <- mex(wavesurge,which=1,mqu=0.7) 
  wavesurge.boot <- bootmex(wavesurge.fit,R=50,trace=51)
  par(mfrow=c(3,2),pty="m")
  check1 <- plot(wavesurge.boot,main="Marginal parameters\nWave surge data of Coles 2001")
  check2 <- plot(wavesurge.boot,plots="dep",main="Dependence parameters\nWave surge data of Coles 2001\nLaplace margins")
  
  # 5-d air pollution data
  
  Qu <- 0.7
  mqus <- c(.9, .7, .7, .85, .7)
  mquw <- 0.7
  smarmex.O3   <- mex(summer, mqu=mqus, which = 1, dqu = Qu, penalty="none",margins="gumbel",constrain=FALSE)
  wmarmex.O3   <- mex(winter, mqu=mquw, which = 1, dqu = Qu, penalty="none",margins="gumbel",constrain=FALSE)
  smarmex.NO2  <- mex(summer, mqu=mqus, which = 2, dqu = Qu, penalty="none",margins="gumbel",constrain=FALSE)
  wmarmex.NO2  <- mex(winter, mqu=mquw, which = 2, dqu = Qu, penalty="none",margins="gumbel",constrain=FALSE)
  smarmex.NO   <- mex(summer, mqu=mqus, which = 3, dqu = Qu, penalty="none",margins="gumbel",constrain=FALSE)
  wmarmex.NO   <- mex(winter, mqu=mquw, which = 3, dqu = Qu, penalty="none",margins="gumbel",constrain=FALSE)
  smarmex.SO2  <- mex(summer, mqu=mqus, which = 4, dqu = Qu, penalty="none",margins="gumbel",constrain=FALSE)
  wmarmex.SO2  <- mex(winter, mqu=mquw, which = 4, dqu = Qu, penalty="none",margins="gumbel",constrain=FALSE)
  smarmex.PM10 <- mex(summer, mqu=mqus, which = 5, dqu = Qu, penalty="none",margins="gumbel",constrain=FALSE)
  wmarmex.PM10 <- mex(winter, mqu=mquw, which = 5, dqu = Qu, penalty="none",margins="gumbel",constrain=FALSE)
  
  Qu <- 0.7
  R <- 50
  
  Sboot.O3 <- bootmex(smarmex.O3, R=R,trace=R+1)
  Wboot.O3 <- bootmex(wmarmex.O3, R=R,trace=R+1)
  Sboot.NO2 <- bootmex(smarmex.NO2, R=R,trace=R+1)
  Wboot.NO2 <- bootmex(wmarmex.NO2, R=R,trace=R+1)
  Sboot.NO <- bootmex(smarmex.NO, R=R,trace=R+1)
  Wboot.NO <- bootmex(wmarmex.NO, R=R,trace=R+1)
  Sboot.SO2 <- bootmex(smarmex.SO2, R=R,trace=R+1)
  Wboot.SO2 <- bootmex(wmarmex.SO2, R=R,trace=R+1)
  Sboot.PM10 <- bootmex(smarmex.PM10, R=R,trace=R+1)
  Wboot.PM10 <- bootmex(wmarmex.PM10, R=R,trace=R+1)
  
  par(mfrow=c(4,2))
  check3 <- plot(Sboot.O3,plots="dep",main="Summer air pollution data\nFig.5 Heffernan & Tawn 2004",xlim=c(0,1),ylim=c(-1,1))
  plot(Wboot.O3,plots="dep",main="Winter air pollution data\nFig.5 Heffernan & Tawn 2004",xlim=c(0,1),ylim=c(-1,1))
  plot(Sboot.NO2,plots="dep",main="Summer air pollution data\nFig.5 Heffernan & Tawn 2004",xlim=c(0,1),ylim=c(-1,1))
  plot(Wboot.NO2,plots="dep",main="Winter air pollution data\nFig.5 Heffernan & Tawn 2004",xlim=c(0,1),ylim=c(-1,1))
  plot(Sboot.NO,plots="dep",main="Summer air pollution data\nFig.5 Heffernan & Tawn 2004",xlim=c(0,1),ylim=c(-1,1))
  plot(Wboot.NO,plots="dep",main="Winter air pollution data\nFig.5 Heffernan & Tawn 2004",xlim=c(0,1),ylim=c(-1,1))
  plot(Sboot.SO2,plots="dep",main="Summer air pollution data\nFig.5 Heffernan & Tawn 2004",xlim=c(0,1),ylim=c(-1,1))
  plot(Wboot.SO2,plots="dep",main="Winter air pollution data\nFig.5 Heffernan & Tawn 2004",xlim=c(0,1),ylim=c(-1,1))
  plot(Sboot.PM10,plots="dep",main="Summer air pollution data\nFig.5 Heffernan & Tawn 2004",xlim=c(0,1),ylim=c(-1,1))
  plot(Wboot.PM10,plots="dep",main="Winter air pollution data\nFig.5 Heffernan & Tawn 2004",xlim=c(0,1),ylim=c(-1,1))
  
  checkEquals(check1,NULL,msg="plot.bootmex successful execution of plotting code 2-d data")
  checkEquals(check2,NULL,msg="plot.bootmex successful execution of plotting code 2-d data")
  checkEquals(check3,NULL,msg="plot.bootmex successful execution of plotting code 3-d data")
} 

test.plot.evmOpt <- function(){
  par(mfrow=c(2,2))
  mod <- evm(rain, th=30, penalty="none")
  res <- plot(mod,main=paste(rep("Figure 4.5 of Coles (2001)",4),
                             c("\nProbability plot","\nQuantile Plot","\nReturn Level Plot\n(SCALE IS DAYS NOT YEARS)","\nDensity Plot")), 
              RetPeriodRange=c(3.65,365*10000))
  checkEquals(res,NULL,msg="plot.evmOpt: GPD successful execution")
  
  # check for very short tailed data
  set.seed(6)
  temp <- rgpd(1000,sigma=1,xi=-0.45)
  fit <- evm(temp,th=0)
  res <- plot(fit,main=c("GPD: PP","GPD: QQ","GPD: RL","GPD: Hist, Short tailed data"))
  checkEquals(res,NULL,msg="plot.evmOpt: GPD successful execution, short tailed data")
  
  # check for covariates in the model
  # GPD
  n <- 1000
  sig <- 2
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,sig,X[,2])
  X$Y <- Y
  fit <- evm(Y,data=X,xi=~b,th=0)
  res <- plot(fit)
  checkEquals(res,NULL,msg="plot.evmOpt: GPD with covariates successful execution")
  
  #GEV 
  # no covariates
  n <- 1000
  Y <- rgev(n,1,1,-.1)
  fit <- evm(Y,family=gev)
  par(mfrow=c(2,2))
  res <- plot(fit,main="GEV no covariates, neg xi")
  checkEquals(res,NULL,msg="plot.evmOpt: GEV no covariates, neg xi successful execution")
  
  Y <- rgev(n,1,1,.2)
  fit <- evm(Y,family=gev)
  res <- plot(fit,main="GEV no covariates, pos xi")
  checkEquals(res,NULL,msg="plot.evmOpt: GEV no covariates, pos xi successful execution")
  
  #GEV with covariates
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3),C= runif(n))
  Y <- rgev(n,X[,3],sig,X[,2])
  X$Y <- Y
  fit <- evm(Y,data=X,xi=~b,mu=~C,family=gev)
  par(mfrow=c(2,2))
  res <- plot(fit,main=rep("GEV with covariates",4))
  checkEquals(res,NULL,msg="plot.evmOpt: GEV with covariates successful execution")
  
  fit <- evm(Y,data=X,xi=~b,family=gev)
  par(mfrow=c(2,2))
  res <- plot(fit,main=rep("GEV with one covariate",4))
  checkEquals(res,NULL,msg="plot.evmOpt: GEV with one covariate successful execution")
}

test.plot.mex <- function(){
  
  smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
  wmarmod <- migpd(winter, mqu=.7,  penalty="none")
  
  mySdep1.lap <- mexDependence(smarmod,which=1, dqu=0.7)
  myWdep1.lap <- mexDependence(wmarmod,which=1, dqu=0.75)
  mySdep2.lap <- mexDependence(smarmod,which=2, dqu=0.8)
  myWdep2.lap <- mexDependence(wmarmod,which=2, dqu=0.85)
  mySdep3.lap <- mexDependence(smarmod,which=3, dqu=0.9)
  myWdep3.lap <- mexDependence(wmarmod,which=3, dqu=0.7)
  mySdep4.lap <- mexDependence(smarmod,which=4, dqu=0.75)
  myWdep4.lap <- mexDependence(wmarmod,which=4, dqu=0.8)
  mySdep5.lap <- mexDependence(smarmod,which=5, dqu=0.85)
  myWdep5.lap <- mexDependence(wmarmod,which=5, dqu=0.9)
  
  mySdep1.gum <- mexDependence(smarmod,which=1, dqu=0.7, margins = "gumbel",constrain=FALSE)
  myWdep1.gum <- mexDependence(wmarmod,which=1, dqu=0.75, margins = "gumbel",constrain=FALSE)
  mySdep2.gum <- mexDependence(smarmod,which=2, dqu=0.8, margins = "gumbel",constrain=FALSE)
  myWdep2.gum <- mexDependence(wmarmod,which=2, dqu=0.85, margins = "gumbel",constrain=FALSE)
  mySdep3.gum <- mexDependence(smarmod,which=3, dqu=0.9, margins = "gumbel",constrain=FALSE)
  myWdep3.gum <- mexDependence(wmarmod,which=3, dqu=0.7, margins = "gumbel",constrain=FALSE)
  mySdep4.gum <- mexDependence(smarmod,which=4, dqu=0.75, margins = "gumbel",constrain=FALSE)
  myWdep4.gum <- mexDependence(wmarmod,which=4, dqu=0.8, margins = "gumbel",constrain=FALSE)
  mySdep5.gum <- mexDependence(smarmod,which=5, dqu=0.85, margins = "gumbel",constrain=FALSE)
  myWdep5.gum <- mexDependence(wmarmod,which=5, dqu=0.9, margins = "gumbel",constrain=FALSE)
  
  
  # check plots produced for variety of thresholds and parameter combinations
  par(mfcol=c(3,4),pty="m")
  plot(mySdep1.gum,main=paste("Summer\nfitting quantile =",mySdep1.gum$dependence$dqu,"\nGumbel margins"))
  plot(myWdep1.gum,main=paste("Winter\nfitting quantile =",myWdep1.gum$dependence$dqu,"\nGumbel margins"))
  plot(mySdep2.gum,main=paste("Summer\nfitting quantile =",mySdep2.gum$dependence$dqu,"\nGumbel margins"))
  plot(myWdep2.gum,main=paste("Winter\nfitting quantile =",myWdep2.gum$dependence$dqu,"\nGumbel margins"))
  plot(mySdep3.gum,main=paste("Summer\nfitting quantile =",mySdep3.gum$dependence$dqu,"\nGumbel margins"))
  plot(myWdep3.gum,main=paste("Winter\nfitting quantile =",myWdep3.gum$dependence$dqu,"\nGumbel margins"))
  plot(mySdep4.gum,main=paste("Summer\nfitting quantile =",mySdep4.gum$dependence$dqu,"\nGumbel margins"))
  plot(myWdep4.gum,main=paste("Winter\nfitting quantile =",myWdep4.gum$dependence$dqu,"\nGumbel margins"))
  plot(mySdep5.gum,main=paste("Summer\nfitting quantile =",mySdep5.gum$dependence$dqu,"\nGumbel margins"))
  plot(myWdep5.gum,main=paste("Winter\nfitting quantile =",myWdep5.gum$dependence$dqu,"\nGumbel margins"))
  
  plot(mySdep1.lap,main=paste("Summer\nfitting quantile =",mySdep1.lap$dependence$dqu,"\nLaplace margins"))
  plot(myWdep1.lap,main=paste("Winter\nfitting quantile =",myWdep1.lap$dependence$dqu,"\nLaplace margins"))
  plot(mySdep2.lap,main=paste("Summer\nfitting quantile =",mySdep2.lap$dependence$dqu,"\nLaplace margins"))
  plot(myWdep2.lap,main=paste("Winter\nfitting quantile =",myWdep2.lap$dependence$dqu,"\nLaplace margins"))
  plot(mySdep3.lap,main=paste("Summer\nfitting quantile =",mySdep3.lap$dependence$dqu,"\nLaplace margins"))
  plot(myWdep3.lap,main=paste("Winter\nfitting quantile =",myWdep3.lap$dependence$dqu,"\nLaplace margins"))
  plot(mySdep4.lap,main=paste("Summer\nfitting quantile =",mySdep4.lap$dependence$dqu,"\nLaplace margins"))
  plot(myWdep4.lap,main=paste("Winter\nfitting quantile =",myWdep4.lap$dependence$dqu,"\nLaplace margins"))
  plot(mySdep5.lap,main=paste("Summer\nfitting quantile =",mySdep5.lap$dependence$dqu,"\nLaplace margins"))
  plot(myWdep5.lap,main=paste("Winter\nfitting quantile =",myWdep5.lap$dependence$dqu,"\nLaplace margins"))
  
  op <- options()
  options(show.error.messages=FALSE)
  checkException(plot.mex(smarmod),msg="plot.mex: exception handle")
  checkException(plot.mex(TRUE),msg="plot.mex: exception handle")
  options(op)
  
  # check execution for 2-d data
  wavesurge.fit <- migpd(wavesurge,mqu=0.8)
  wavesurge.mex <- mexDependence(wavesurge.fit,dqu=0.8,which=2,margins="gumbel",constrain=FALSE)
  par(mfrow=c(2,2))
  res <- plot(wavesurge.mex,main="Wave surge data")
  checkEquals(res,NULL,msg = "plot.mex: successful execution")
}

test.plot.predict.mex <- function(){
  # check reproduce Figure 6 in Heffernan and Tawn
  w <- mex(winter,mqu=0.7,penalty="none", which="NO", dqu=.7, margins="gumbel", constrain=FALSE)
  noMod <- bootmex(w,trace=101)
  noPred <- predict(noMod,trace=101)
  par(mfcol=c(2,2))
  res <- plot(noPred,main="Fig. 6 Heffernan and Tawn (2004)")
  checkEquals(res, NULL, msg="plot.predict.mex: correct execution")
  
  # check for 2-d data
  R <- 20
  nsim <- 100
  wavesurge.mex <- mex(wavesurge,mqu=.7,which=1, margins="laplace")
  wavesurge.boot <- bootmex(wavesurge.mex,R=R,trace=R+1)
  wavesurge.pred <- predict(wavesurge.boot,nsim=nsim,trace=R+1)
  par(mfrow=c(1,1))
  res <- plot(wavesurge.pred)
  checkEquals(res, NULL, msg="plot.predict.mex: correct execution")
}

test.plot.migpd <- function(){
  par(mfrow=c(2,2))
  mod <- migpd(winter, mqu=.7, penalty = "none")
  res <- plot(mod)
  checkEquals(res,NULL,msg="plot.migpd: successful execution")
}

test.plot.predict.evm <- function(){
  # testing all of: plot.lp.evm* and plot.rl.evm* where * is opt, sim and boot
  
  # first with no covariates
  n <- 100
  for(Family in list(gpd,gev)){
    set.seed(20130513)
    pst <- function(msg) texmexPst(msg,Family=Family)
    
    u    <- switch(Family$name,GPD=14,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    fit.opt <- evm(data,th=u,family=Family)
    fit.sim <- evm(data,th=u,method="sim",trace=100000,family=Family)
    fit.boot <- evmBoot(fit.opt,R=20,trace=30)
    
    M <- seq(20,1000,len=20)
    
    p.opt <- predict(fit.opt,M=M,ci=TRUE)
    p.sim <- predict(fit.sim,M=M,ci=TRUE)
    p.boot <- predict(fit.boot,M=M,ci=TRUE)
    
    par(mfrow=c(3,3))
    plot(p.opt,main=paste(Family$name,"\nMLE"))
    plot(p.sim,main=paste(Family$name,"\nMCMC"))
    plot(p.boot,main=paste(Family$name,"\nBootstrap"))
    
    p.lp.opt <- predict(fit.opt,type="lp",ci=TRUE)
    p.lp.sim <- predict(fit.sim,type="lp",ci=TRUE)
    p.lp.boot <- predict(fit.boot,type="lp",ci=TRUE)
    
    checkException(plot(p.lp.opt),silent=TRUE,msg=pst("plot.lp.evmOpt: fail if no covariates"))
    checkException(plot(p.lp.sim),silent=TRUE,msg=pst("plot.lp.evmSim: fail if no covariates"))
    checkException(plot(p.lp.boot),silent=TRUE,msg=pst("plot.lp.evmBoot: fail if no covariates"))
    
    # now with covariates
    
    n <- 1000
    M <- 1000
    
    mu <- 1
    
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=cbind(X[,1],X[,2]),GEV=cbind(mu,X[,1],X[,2]))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    start <- switch(Family$name,GPD=c(0,1,0,1),GEV=c(mu,0,1,0,1))
    
    fit.opt <- evm(Y,data=X,phi=~a,xi=~b, th=th,family=Family,start=start)
    fit.sim <- evm(Y,data=X,phi=~a,xi=~b, th=th,family=Family,method="sim",trace=100000,start=start)
    o <- options(warn=-1)
    fit.boot <- evmBoot(fit.opt,R=20,trace=30)
    options(o)
    
    nx <- 3
    M <- seq(5,1000,len=20)
    newX <- data.frame(a=rnorm(nx),b=runif(nx,-0.1,0.1))
    
    p.opt <- predict(fit.opt,M=M,newdata=newX,ci=TRUE)
    p.sim <- predict(fit.sim,M=M,newdata=newX,ci=TRUE)
    p.boot <- predict(fit.boot,M=M,newdata=newX,ci=TRUE)
    
    p.lp.opt <- predict(fit.opt,type="lp",newdata=newX,ci=TRUE)
    p.lp.sim <- predict(fit.sim,type="lp",newdata=newX,ci=TRUE)
    p.lp.boot <- predict(fit.boot,type="lp",newdata=newX,ci=TRUE)
    
    par(mfrow=c(3,3))
    plot(p.opt,sameAxes=FALSE,main=paste(Family$name,"MLE\ndifferent axes"))
    plot(p.sim,sameAxes=FALSE,main=paste(Family$name,"MCMC\ndifferent axes"))
    plot(p.boot,sameAxes=FALSE,main=paste(Family$name,"Bootstrap\ndifferent axes"))
    
    plot(p.opt,sameAxes=TRUE,main=paste(Family$name,"MLE\nsame axes"))
    plot(p.sim,sameAxes=TRUE,main=paste(Family$name,"MCMC\nsame axes"))
    plot(p.boot,sameAxes=TRUE,main=paste(Family$name,"Bootstrap\nsame axes"))
    
    par(mfrow=c(4,4))
    plot(p.lp.opt,main=paste(Family$name,"MLE"))
    plot(p.lp.sim,main=paste(Family$name,"MCMC\nmean"),type="mean")
    plot(p.lp.sim,main=paste(Family$name,"MCMC\nmedian"),type="median")
    plot(p.lp.boot,main=paste(Family$name,"Bootstrap"))
    
    # single covariate only:
    
    param <- switch(Family$name,GPD=cbind(X[1,1],X[,2]),GEV=cbind(mu,X[1,1],X[,2]))
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit.opt <- evm(Y,data=X,xi=~b,th=th,family=Family)
    fit.sim <- evm(Y,data=X,xi=~b,th=th,family=Family,method="sim",trace=100000)
    o <- options(warn=-1)
    fit.boot <- evmBoot(fit.opt,R=20,trace=30)
    options(o)
    
    p.opt <- predict(fit.opt,M=M,newdata=newX,ci=TRUE)
    p.sim <- predict(fit.sim,M=M,newdata=newX,ci=TRUE)
    p.boot <- predict(fit.boot,M=M,newdata=newX,ci=TRUE)
    
    p.lp.opt <- predict(fit.opt,type="lp",newdata=newX,ci=TRUE)
    p.lp.sim <- predict(fit.sim,type="lp",newdata=newX,ci=TRUE)
    p.lp.boot <- predict(fit.boot,type="lp",newdata=newX,ci=TRUE)
    
    par(mfrow=c(3,3))
    plot(p.opt,sameAxes=FALSE,main=paste(Family$name,"MLE"))
    plot(p.sim,sameAxes=FALSE,main=paste(Family$name,"MCMC"))
    plot(p.boot,sameAxes=FALSE,main=paste(Family$name,"Bootstrap"))
    
    par(mfrow=c(3,1))
    plot(p.lp.opt,main=paste(Family$name,"MLE"),polycol="cyan")
    plot(p.lp.sim,main=paste(Family$name,"MCMC"),polycol="cyan")
    plot(p.lp.boot,main=paste(Family$name,"Bootstrap"),polycol="cyan")
  }
}

test.plotrl.evm <- function(){
  # no covariates
  
  rain.fit <- evm(rain,th=30)
  par(mfrow=c(1,1))
  plotrl.evmOpt(rain.fit,RetPeriodRange=c(1,2000),main="Coles (2001) figure 4.5\nGPD Return Level Plot")
  
  sealevel.fit <- evm(portpirie$SeaLevel,family=gev)
  plotrl.evmOpt(sealevel.fit,main="Coles (2001), Figure 3.5\nGEV Return Level Plot")
  
  # with covariates
  
  for(Family in list(gpd,gev)){
    pst <- function(msg) texmexPst(msg,Family=Family)
    set.seed(20130513)
    n <- 100    
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=X,GEV=cbind(5,X))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,phi=~a,xi=~b,th=th,family=Family)
    
    checkException(plotrl.evmOpt(fit),silent=TRUE,msg=pst("plotrl.evmOpt : failure for model with covariates"))
  }
}

################################################################################
## test.predict.evm()

test.predict.evmOpt <- function(){
  
  for(Family in list(gpd,gev)){
    set.seed(20130513)
    pst <- function(msg) texmexPst(msg,Family=Family)
    
    # no covariates
    
    u    <- switch(Family$name,GPD=14,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    r.fit <- evm(data,th=u,family=Family)
    co <- coef(r.fit)
    
    if(Family$name == "GPD")checkEqualsNumeric(target=u,current = predict(r.fit,M=1/r.fit$rate)[[1]],msg=pst("predict.evmOpt: GPD retrieve threshold"))
    
    checkEquals(target=predict(r.fit), current=rl(r.fit),msg=pst("predict.evmOpt: predict with type=rl gives same as direct call to rl with default arguments"))
    checkEquals(target=predict(r.fit,type="lp"), current=linearPredictors(r.fit),msg=pst("predict.evmOpt: predict with type=lp gives same as direct call to linearpredictors with default arguments"))
    
    r.fit$rate <- 1
    prob <- c(0.5,0.9,0.95,0.99,0.999)
    for(p in prob){
      checkEqualsNumeric(target = Family$quant(p,t(co),r.fit),
                         current = unlist(predict(r.fit,M=1/(1-p))),msg=pst("predict.evmOpt: est ret levels no covariates"))
    }
    
    # with covariates
    
    n <- 1000
    M <- 1000
    
    # boundary cases with single covariates in only one parameter
    mu <- 1
    xi <- 0.05
    
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=cbind(X[,1],xi),GEV=cbind(mu,X[,1],xi))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,phi=~a,th=th,family=Family)
    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]
    pred <- predict(fit,M=M)
    
    checkEqualsNumeric(target=X$a,current=pred[[1]][,-1],msg=pst("predict.evmOpt: ret level correct reporting of covariates for single cov in phi only"))
    checkEqualsNumeric(target=Family$quant(1-1/M,AllCo,fit),
                       current = pred[[1]][,1],msg=pst("predict.evmOpt: ret level estimation with covariates in phi only"))
    
    mu <- 1
    sig <- 2
    
    param <- switch(Family$name,GPD=cbind(log(sig),X[,2]),GEV=cbind(mu,log(sig),X[,2]))
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,xi=~b,th=th,family=Family)
    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]
    pred <- predict(fit,M=M)
    
    checkEqualsNumeric(target=X$b,current=pred[[1]][,-1],msg=pst("predict.evmOpt: ret level correct reporting of covariates for single cov in xi only"))
    checkEqualsNumeric(target=Family$quant(1-1/M,AllCo,fit),
                       current = pred[[1]][,1],msg=pst("predict.evmOpt: ret level estimation with covariates in xi only"))
    
    # covariates in all parameters
    
    param <- switch(Family$name,GPD=cbind(X[,1],X[,2]),GEV=cbind(X[,1],X[,1],X[,2]))
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- switch(Family$name,
                  GPD = evm(Y, data=X,        phi=~a, xi=~b, th=th,family=Family),
                  GEV = evm(Y, data=X, mu=~a, phi=~a, xi=~b, th=th,family=Family))  
    
    AllCo <- predict(fit,type="lp")$link[,1:length(Family$param)]
    
    checkEqualsNumeric(target=Family$quant(1-1/M,AllCo,fit),
                       current = predict(fit,M=M)[[1]][,1],msg=pst("predict.evmOpt: ret level estimation with covariates in all parameters"))
    
    # check multiple M
    M <- c(10,50,100,500,1000)
    
    target <- sapply(M,function(m,AllCo,fit) Family$quant(1-1/m,AllCo,fit),AllCo=AllCo,fit=fit)
    current <- predict(fit,M=M)
    
    for(i in 1:length(M)){
      checkEqualsNumeric(target[,i],current[[i]][,1],msg=pst("predict.evmOpt: ret level estimation multiple M"))
    }
    
    # new data
    nx <- 20
    M <- 1000
    newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
    AllCoNew <- predict(fit,type="lp",newdata=newX)$link[,1:length(Family$param)]
    
    checkEqualsNumeric(target=Family$quant(1-1/M,AllCoNew,fit),current=predict(fit,M=M,newdata=newX)[[1]][,1],msg=pst("predict.evmOpt: ret level ests with new data"))
    
    checkEqualsNumeric(dim(AllCoNew) + c(0,3),dim(predict(fit,ci=TRUE,newdata=newX)[[1]]), msg=pst("predict.evmOpt: dimension of return object for ci calc"))
    checkEqualsNumeric(dim(AllCoNew) + c(0,2),dim(predict(fit,se=TRUE,newdata=newX)[[1]]), msg=pst("predict.evmOpt: dimension of return object for se calc"))
    checkEqualsNumeric(dim(AllCoNew) + c(0,4),dim(predict(fit,se=TRUE,ci=TRUE,newdata=newX)[[1]]), msg=pst("predict.evmOpt: dimension of return object for se and ci calc"))
    
    Labels <- switch(Family$name,GPD=c("a","b"),GEV=c("a","a","b"))
    
    checkEquals(c("RL","2.5%","97.5%","se.fit",Labels), colnames(predict(fit,se=TRUE,ci=TRUE)[[1]]), msg=pst("predict.evmOpt: colnames of return obejct for se and ci calc, default alpha"))
    checkEquals(c("RL","5%","95%","se.fit",Labels), colnames(predict(fit,se=TRUE,ci=TRUE,alpha=0.1)[[1]]), msg=pst("predict.evmOpt: colnames of return obejct for se and ci calc, alpha=0.1"))
    
    # alpha
    
    alpha <- c(0.01,0.05,0.1,0.2,0.5,0.9,0.99)
    z <- matrix(qnorm(c(alpha/2,1-alpha/2)),ncol=2)
    
    for(a in 1:length(alpha)){
      p <- predict(fit,alpha=alpha[a],ci=TRUE,se=TRUE)[[1]]
      checkEquals(current = colnames(p)[2:3],target = c(paste(100*alpha[a]/2,"%",sep=""),paste(100*(1-alpha[a]/2),"%",sep="")),msg=pst("predict.evmOpt: labelling of confidence intervals"))
      checkEqualsNumeric(target = p[,1] + z[a,1]*p[,4],current = p[,2], msg=pst("predict.evmOpt: ret level Conf Interval calc for different alpha"))
      checkEqualsNumeric(target = p[,1] + z[a,2]*p[,4],current = p[,3], msg=pst("predict.evmOpt: ret level Conf Interval calc for different alpha"))
    }
    
    
    # linear predictors
    
    checkEqualsNumeric(target = c(nx,switch(Family$name,GPD=6,GEV=9)), dim(predict(fit,newdata=newX,se=TRUE,type="lp")$link), msg=pst("predict.evmOpt: dimension of return object, linear predictor, se calc"))
    checkEqualsNumeric(target = c(nx,switch(Family$name,GPD=8,GEV=12)),dim(predict(fit,newdata=newX,ci=TRUE,type="lp")$link), msg=pst("predict.evmOpt: dimension of return object, linear predictor, ci calc"))
    
    nameCI <- switch(Family$name,GPD = c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi","a","b"),GEV = c("mu","phi", "xi", "mu.lo","mu.hi","phi.lo", "phi.hi", "xi.lo", "xi.hi","a","a","b"))
    nameSE <- switch(Family$name,GPD = c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi","phi.se", "xi.se","a","b"),GEV = c("mu","phi", "xi", "mu.lo","mu.hi","phi.lo", "phi.hi", "xi.lo", "xi.hi","mu.se","phi.se", "xi.se","a","a","b"))
    checkEquals(target = nameCI, current = colnames(predict(fit,newdata=newX,ci=TRUE,type="lp")$link),msg=pst("predict.evmOpt: colnames for linear predictor return object"))
    checkEquals(target = nameSE, current = colnames(predict(fit,newdata=newX,ci=TRUE,se=TRUE,type="lp")$link),msg=pst("predict.evmOpt: colnames for linear predictor return object"))
    
    # unique
    
    newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1)) # checks for duplicates in one and both covariates.
    U <- !duplicated(newX)
    checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")$link,
                       target = predict(fit,newdata=newX,unique.=FALSE,type="lp")$link[U,], msg=pst("predict.evmOpt: functioning of argument unique, for linear predictor"))
    checkEqualsNumeric(current = predict(fit,newdata=newX)[[1]],
                       target =  predict(fit,newdata=newX,unique.=FALSE)[[1]][U,], msg=pst("predict.evmOpt: functioning of argument unique, for return levels"))
    
    # check standard errors - this takes a while since using bootstrap
    
    M <- c(10,100,500,1000,2000)
    newX <- data.frame("a"=rep(c(1,-1,2,-2),2),"b"=c(rep(0.1,4),rep(-0.1,4)))
    fit.p <- predict(fit, newdata=newX,se=TRUE,M=M)
    fit.seest <- unlist(lapply(fit.p,function(x) x[,2]))
    
    o <- options(warn=-1)
    fit.b <- evmBoot(fit,R=1000, trace=1100)
    options(o)
    fit.bp <- predict(fit.b,newdata=newX,all=TRUE,M=M)
    fit.seb <- lapply(fit.bp,function(X) apply(X,2,sd))
    fit.seboot <- unlist(fit.seb)
    
    checkTrue(all(abs((fit.seboot -  fit.seest) / fit.seest) < 0.3),msg=pst("predict.evmOpt: return level standard error estimate compared with bootstrap standard errors"))
  }  
}

################################################################################
## test.predict.evmSim()

test.predict.evmSim <- function(){
  
  for(Family in list(gpd,gev)){
    
    pst <- function(msg) texmexPst(msg,Family=Family)
    set.seed(20130513)
    # no covariates
    
    u    <- switch(Family$name,GPD=14,GEV=-Inf)
    data <- switch(Family$name,GPD=rain,GEV=portpirie$SeaLevel)
    
    r.fit <- evm(data,th=u,family=Family,method="sim",trace=50000)
    co <- coef(r.fit)
    
    if(Family$name == "GPD"){
      checkEqualsNumeric(target=u,current=predict(r.fit,M=1/r.fit$map$rate)[[1]], msg=pst("predict.evmSim: retrieve threshold"))
    }
    
    r.fit$map$rate <- 1
    p <- c(0.5,0.9,0.95,0.99,0.999)
    checkEqualsNumeric(target = Family$quant(p,t(co),r.fit$map), tolerance=0.01,
                       current = unlist(predict(r.fit,M=1/(1-p))),msg=pst("predict.evmSim: ret level estimation"))
    
    checkEquals(target=predict(r.fit,M=1/(1-p)), current=rl(r.fit,M=1/(1-p)),msg=pst("predict.evmSim: predict with type=rl gives same as direct call to rl with default arguments"))
    checkEquals(target=predict(r.fit,type="lp")$link, current=linearPredictors(r.fit)$link,msg=pst("predict.evmSim: predict with type=rl gives same as direct call to rl with default arguments"))
    
    # with covariates
    
    n <- 1000
    M <- 1000
    mu <- 1
    
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=X,GEV=cbind(mu,X))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    start <- switch(Family$name,GPD=c(0,1,0,1),GEV=c(1,0,1,0,1))
    fit <- evm(Y,data=X,phi=~a,xi=~b,th=th,method="sim",trace=50000,family=Family,start=start)
    
    AllCo <- predict(fit,type="lp",all=TRUE)$link
    PostMeanRL <- function(AllCo,M){
      AllQuant <- lapply(AllCo,function(X)Family$quant(1-1/M,X[,1:length(Family$param)],fit$map))
      sapply(AllQuant,mean)
    }
    
    checkEqualsNumeric(target=PostMeanRL(AllCo,M),current = predict(fit,M=M)[[1]][,1],msg=pst("predict.evmSim: ret level estimation with covariates in all parameters"))
    # multiple M
    
    M <- c(10,50,100,500,1000)
    
    target <- lapply(M,function(m) PostMeanRL(AllCo,m))
    current <- predict(fit,M=M)
    
    for(i in 1:length(M)){
      checkEqualsNumeric(target[[i]],current[[i]][,1],tolerance=0.02,msg=pst("predict.evmSim: ret level estimation multiple M"))
    }
    
    # new data
    nx <- 20
    M <- 1000
    newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
    
    AllCoNew <- predict(fit,type="lp",newdata=newX,all=TRUE)$link
    
    checkEqualsNumeric(target=PostMeanRL(AllCoNew,M),current=predict(fit,M=M,newdata=newX)[[1]][,1],msg=pst("predict.evmSim: ret level ests with new data"))
    checkEqualsNumeric(target = as.matrix(newX), current = predict(fit,M=M,newdata=newX)[[1]][,2:3],msg=pst("predict.evmSim: ret level estimation with new data, covariates added correctly to output"))
    
    
    p <- predict(fit,all=TRUE,newdata=newX)
    alpha <- c(0.05,0.1)
    for(a in alpha){
      l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
      l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
      m.L <- unlist(l.L)
      m.U <- unlist(l.U)
      r <- predict(fit,newdata=newX,ci=TRUE,alpha=a)
      checkEqualsNumeric(target = m.L,current = r[[1]][,3],msg=pst("predict.evmSim : lower conf ints for ret levels with new data"))
      checkEqualsNumeric(target = m.U,current = r[[1]][,4],msg=pst("predict.evmSim : upper conf ints for ret levels with new data"))
    }
    
    # check linear predictors
    p <- predict(fit,type="lp",all=TRUE,newdata=newX)$link
    l <- lapply(p,function(l) apply(l,2,mean))
    m <- matrix(unlist(l),ncol=length(l[[1]]),byrow=TRUE)
    r <- predict(fit,type="lp",newdata=newX)$link
    checkEqualsNumeric(target = m,current = r,msg=pst("predict.evmSim : linear predictors of parameters with new data"))
    Offset <- switch(Family$name,GEV=1,GPD=0)
    checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")$link[,Offset +1:2],
                       target = cbind(apply(cbind(rep(1,nx),newX[,1]) %*% t(fit$param[,Offset + 1:2]),1,mean),
                                      apply(cbind(rep(1,nx),newX[,2]) %*% t(fit$param[,Offset + 3:4]),1,mean)),
                       msg = pst("predict.evmSim: linear predictor estimates"))
    
    alpha <- c(0.05,0.1)
    for(a in alpha){
      l.L <- lapply(p,function(l) apply(l,2,quantile,prob=a/2))
      l.U <- lapply(p,function(l) apply(l,2,quantile,prob=1-a/2))
      m.L <- matrix(unlist(l.L),ncol=length(l[[1]]),byrow=TRUE)
      m.U <- matrix(unlist(l.U),ncol=length(l[[1]]),byrow=TRUE)
      r <- predict(fit,type="lp",newdata=newX,ci=TRUE,alpha=a)$link
      npar <- length(Family$param)
      checkEqualsNumeric(target = m.L[,1:npar],current = r[,(1:(4*npar))[rep(c(F,F,T,F),npar)]],msg=pst("predict.evmSim : lower conf ints for linear predictors of parameters with new data"))
      checkEqualsNumeric(target = m.U[,1:npar],current = r[,(1:(4*npar))[rep(c(F,F,F,T),npar)]],msg=pst("predict.evmSim : upper conf ints for linear predictors of parameters with new data"))
    }
    
    # structure of output
    
    checkEqualsNumeric(target = c(n,6), current = dim(predict(fit,ci=TRUE)[[1]]), msg=pst("predict.evmSim: dimension of output with ci calculation"))
    o <- options(warn=-1) # since se=TRUE gives a warning
    checkEqualsNumeric(target = c(n,6), current = dim(predict(fit,se=TRUE,ci=TRUE)[[1]]), msg=pst("predict.evmSim: dimension of output with ci and se calculation"))
    options(o)
    
    checkEquals(target = c("Mean", "50%","2.5%","97.5%","a","b"), colnames(predict(fit,ci=TRUE)[[1]]), msg=pst("predict.evmSim: colnames of ret level ests with CI estimation"))
    checkEquals(target = c("Mean", "50%","5%","95%","a","b"), colnames(predict(fit,alpha=0.1,ci=TRUE)[[1]]), msg=pst("predict.evmSim: colnames of ret level ests with CI estimation, alpha=0.1"))
    
    checkEqualsNumeric(target = c(nx,4*npar+2), dim(predict(fit,newdata=newX,ci=TRUE,type="lp")$link), msg=pst("predict.evmSim: dimension of linear predictor return object"))
    
    cnamesGPD <- c("phi: Mean", "phi: 50%", "phi: 2.5%", "phi: 97.5%", "xi: Mean", "xi: 50%", "xi: 2.5%", "xi: 97.5%")# this specific format assumed by plot.rl.evmSim and plot.lp.evmSim
    cnamesGEV <- c("mu: Mean",  "mu: 50%",  "mu: 2.5%",  "mu: 97.5%",  "phi: Mean", "phi: 50%", "phi: 2.5%", "phi: 97.5%", "xi: Mean", "xi: 50%", "xi: 2.5%", "xi: 97.5%")
    cnames<- switch(Family$name,GPD=cnamesGPD,GEV=cnamesGEV)
    
    checkEquals(current = cnames, target = colnames(predict(fit,newdata=newX,ci=TRUE,type="lp")$link)[1:(4*npar)], msg=pst("predict.evmSim: col names of lin predictors with CI calcs"))
    o <- options(warn=-1) # since se=TRUE gives a warning
    checkEquals(current = cnames, target = colnames(predict(fit,newdata=newX,ci=TRUE,se=TRUE,type="lp")$link)[1:(4*npar)], msg=pst("predict.evmSim: col names of lin predictors with CI+SE calcs"))
    options(o)
    
    # unique
    newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1))
    
    checkEqualsNumeric(current = predict(fit,newdata=newX)[[1]], target = unique(predict(fit,newdata=newX,unique.=FALSE)[[1]]),pst("predict.evmSim: unique functioning for ret level ests"))
    checkEqualsNumeric(current = predict(fit,newdata=newX,type="lp")$link[,], target = unique(predict(fit,newdata=newX,unique.=FALSE,type="lp")$link[,]),msg=pst("predict.evmSim: unique functioning for lin pred ests"))
  }
}

################################################################################
## test.predict.evmBoot()

test.predict.evmBoot <- function(){
  # functionality all tested already in test.predict.evm, so just check output of correct format.
  
  n <- 1000
  nx <- 9
  R <- 10
  nm <- 20
  
  for(Family in list(gpd,gev)){
    
    pst <- function(msg) texmexPst(msg,Family=Family)
    set.seed(20130513)
    
    X <- data.frame(a = rnorm(n),b = runif(n,-0.1,0.1))
    param <- switch(Family$name,GPD=X,GEV=cbind(5,X))
    th <- switch(Family$name,GPD=0,GEV=-Inf)
    X$Y <- Family$rng(n,param,list(threshold=th))
    fit <- evm(Y,data=X,phi=~a,xi=~b,th=th,family=Family)
    
    o <- options(warn=-1)
    boot <- evmBoot(fit,R=R,trace=100)
    options(o)
    
    newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
    from <- 10; to <- 500
    M <- seq(from,to,length=nm)
    pred <- predict(boot,newdata=newX,M=M,ci=TRUE)
    
    checkEquals(target=predict(boot), current=rl(boot),msg=pst("predict.evmBoot: predict with type=rl gives same as direct call to rl with default arguments"))
    checkEquals(target=predict(boot,type="lp"), current=linearPredictors(boot),msg=pst("predict.evmBoot: predict with type=lp gives same as direct call to linearPredictors with default arguments"))
    
    checkEqualsNumeric(target=nm,current=length(pred),msg=pst("predict.evmBoot: output length"))
    checkEquals(target=paste("M.",from,sep=""),current=names(pred)[1],msg=pst("predict.evmBoot: names of output"))
    checkEquals(target=paste("M.",to,sep=""),current=names(pred)[nm],msg=pst("predict.evmBoot: names of output"))
    
    cnames <- c( "Mean","50%","2.5%","97.5%",names(X)[1:2])
    checkEquals(target=cnames,current=colnames(pred[[1]]),msg=pst("predict.evmBoot: colnames"))
    
    checkEqualsNumeric(target=c(nx,6),current=dim(pred[[1]]),msg=pst("predict.evmBoot: dimension"))
    for(i in 1:nm){
      checkEqualsNumeric(target=newX[,1],current=pred[[i]][,5],msg=pst("predict.evmBoot: covariates inoutput"))
      checkEqualsNumeric(target=newX[,2],current=pred[[i]][,6],msg=pst("predict.evmBoot: covariates inoutput"))
    }
    
    par(mfrow=n2mfrow(nx))
    plot(pred,sameAxes=FALSE,type="median",main="Bootstrap median rl")
    plot(pred,sameAxes=FALSE,type="mean",main="Bootstrap mean rl")
  }
}

test.predict.mex <- function(){
  # reproduce Table 5 in Heffernan and Tawn 2004
  smarmod <- mex(summer, mqu=c(.9, .7, .7, .85, .7), which="NO", penalty="none", dqu=.7,margins="gumbel",constrain=FALSE)
  wmarmod <- mex(winter, mqu=.7,  penalty="none", which="NO",margins="gumbel",constrain=FALSE)
  set.seed(20111010)
  NOmodWinter <- bootmex(wmarmod,trace=101)
  NOpredWinter <- predict(NOmodWinter, nsim = 500,trace=101) # matches sample size in H+T2004
  
  NOmodSummer <- bootmex(smarmod,trace=101)
  NOpredSummer <- predict(NOmodSummer, nsim = 500,trace=101)
  
  Table5winter <- rbind(c(8.3, 75.4, 569.9, 44.6, 132.3),
                        c(1.2, 4.4, 45.2, 6.7, 8.2))
  Table5summer <- rbind(c(39.6,62.2,213.5,48.5,83.7),
                        c(4.3,4.3,17.5,11.8, 7.9))
  
  dimnames(Table5winter) <- dimnames(Table5summer) <- list(c("E(x)", "SE"),
                                                           c("O3", "NO2", "NO", "SO2", "PM10"))
  
  Table5summer <- Table5summer[, c("NO", "O3", "NO2", "SO2", "PM10")]
  Table5winter <- Table5winter[, c("NO", "O3", "NO2", "SO2", "PM10")]
  
  resSummer <- summary(NOpredSummer)$ans[1:2,]
  resWinter <- summary(NOpredWinter)$ans[1:2,]
  
  pointEstSummer <- apply(NOpredSummer$data$sim,2,mean)
  pointEstWinter <- apply(NOpredWinter$data$sim,2,mean)
  
  tol <- 0.05
  
  checkEqualsNumeric(Table5summer, resSummer,tolerance=tol,msg="predict.mex: Table 5 summer data")
  checkEqualsNumeric(Table5winter, resWinter,tolerance=tol,msg="predict.mex: Table 5 winter data")
  
  checkEqualsNumeric(pointEstSummer, resSummer[1,],tolerance=tol,msg="predict.mex: point est vs boot, summer data")
  checkEqualsNumeric(pointEstWinter, resWinter[1,],tolerance=tol,msg="predict.mex: point est vs boot, winter data")
  
  # check execution for 2-d data
  
  R <- 20
  nsim <- 100
  wavesurge.mex <- mex(wavesurge,mqu=.7,dqu=0.7,margins="laplace",which=1)
  wavesurge.boot <- bootmex(wavesurge.mex,R=R,trace=R+1)
  wavesurge.pred <- predict(wavesurge.boot,nsim=nsim,trace=R+1)
  
  checkEqualsNumeric(length(wavesurge.pred$replicates),R,msg="predict.mex execution for 2-d data")
  checkEqualsNumeric(dim(wavesurge.pred$replicates[[3]]),c(nsim,2))
  checkEquals(names(wavesurge.pred$replicates[[4]]),names(wavesurge),msg="predict.mex execution for 2-d data")
  
  # check predictions Laplace estimation equal to Gumbel for large samples and high threshold
  
  tol <- 0.01
  seeds <- 20:24
  set.seed(20111004)
  n <- 100000
  mqu <- c(0,0.9)
  dqu <- 0.99
  for(i in 1:5){
    x <- rgpd(n=n,sigma=1,xi=0.1)
    y <- 5 + rexp(1,5)*x + rnorm(n,0,x/max(x))
    data <- data.frame(x=x,y=y)
    
    data.gpd <- migpd(data , mqu=mqu, penalty="none")
    lap.mex <- mexDependence(data.gpd,which=1, dqu=dqu,start=c(-0.1,0.1),PlotLikDo=FALSE,v=20)
    gum.mex <- mex(data,mqu=c(0,0.9),which=1, dqu=dqu,margins="gumbel",constrain=FALSE)
    
    set.seed(seeds[i])
    lap.pred <- predict(lap.mex,nsim=10000,trace=R+1)
    set.seed(seeds[i])
    gum.pred <- predict(gum.mex,nsim=10000,trace=R+1)
    
    lap.ans <- summary(lap.pred)$ans
    gum.ans <- summary(gum.pred)$ans
    
    checkEqualsNumeric(lap.ans,gum.ans,tolerance=tol,msg=paste("predict.mex Laplace predictions equal to Gumbel, test replicate",i))
  }
}

test.qgev <- function() {
  ## get the probabilities that we'll use and sort them
  ## into ascending order for safekeeping
  probabilities <- sort(runif(10))
  
  core.sanity.test <- function(xi) {
    base.quantiles <- qgev(probabilities, 0, 1, xi)
    ## check that the values are ascending
    checkTrue(all(diff(base.quantiles) >= 0), "qgev: ascending quantiles")
    ## and check that we're descending correctly for upper tail
    bq2 <- qgev(probabilities, 0, 1, xi, lower.tail=FALSE)
    checkTrue(all(diff(bq2) <= 0), "qgev: descending quantiles")
    ## does lower.tail work
    checkEqualsNumeric(base.quantiles,
                       qgev(1 - probabilities, 0, 1, xi, lower.tail=FALSE),
                       "qgev: lower.tail works correctly")
    ## does log.p work?
    checkEqualsNumeric(base.quantiles,
                       qgev(log(probabilities), 0, 1, xi, log.p=TRUE),
                       "qgev: log.p works")
    ## check shift and scale property
    sigma <- rexp(1)
    mu    <- runif(1, -5, 5)
    shifted <- mu + sigma * base.quantiles
    checkEqualsNumeric(shifted,
                       qgev(probabilities, mu, sigma, xi),
                       "qgev: shift and scale")
  }
  
  lapply(c(0, seq(-5, 5, length.out=10)), core.sanity.test)
  
  ## known values
  checkEqualsNumeric(-log(log(2)),
                     qgev(0.5, 0, 1, 0),
                     "qgev: median match at zero xi")
  xi <- seq(-2, 2, length=10)
  checkEqualsNumeric(qgev(0.5, 0, 1, xi),
                     expm1(-log(log(2))*xi) / xi,
                     "qgev: median match at nonzero xi")
}

test.qgpd <- function(){
  
  set.seed(201110101)
  evd.qgpd <- texmex:::.evd.qgpd
  myTest <- function(sig,xi,thresh,msg){
    myq <- sapply(1:nreps,function(i) qgpd(x[,i], sig[i], xi[i], u=thresh[i]))
    myp <- sapply(1:nreps,function(i) pgpd(myq[,i], sig[i], xi[i], u=thresh[i]))
    eq <- sapply(1:nreps, function(i) evd.qgpd(x[,i], loc=thresh[i], scale=sig[i], shape=xi[i]))
    checkEqualsNumeric(eq,myq,msg=paste(msg,"test using .evd.qgpd"))
    checkEqualsNumeric(x,myp,msg=paste(msg,"test using qgpd"))
  }
  
  #*************************************************************
  # 6.4.0 Test exception for out of range probabilties
  op <- options()
  options(show.error.messages=FALSE)
  checkException(qgpd(1.5,1,0,2),msg="qgpd: exception for out of range prob")
  checkException(qgpd(-1,1,0,2),msg="qgpd: exception for out of range prob")
  options(op)
  
  #*************************************************************
  # 6.4. Test qgpd. Note that .evd.qgpd is NOT vectorized.
  
  nreps <- 100
  nsim <- 1000
  p <- matrix(runif(2*nreps, -1, 1),ncol=2) 
  p[, 1] <- p[, 1] + 1
  thresh <- rep(0,nreps) 
  x <- matrix(runif(nreps*nsim), nrow=nsim)
  
  myTest(sig=p[,1], xi=p[,2],thresh=thresh,msg="qgpd: random xi")
  
  #*************************************************************
  # 6.5. Test qgpd when some or all of xi == 0. Note that .evd.qgpd is NOT vectorized.
  
  p[sample(1:nreps,nreps/2),2] <- 0
  myTest(sig=p[,1], xi = p[,2], thresh=thresh,msg="qgpd: some zero xi")
  p[,2] <-  0
  myTest(sig=p[,1], xi = p[,2], thresh=thresh,msg="qgpd: all zero xi")
  
  #*************************************************************
  # 6.6. Test vectorization of qgpd. Note that .evd.qgpd is NOT vectorized.
  
  sig <- runif(nsim, 0, 2)
  xi <- runif(nsim)
  thresh <- rnorm(nsim)
  
  x <- runif(nsim)
  
  myq <- qgpd(x, sig, xi, thresh)
  eq <- sapply(1:nsim, function(i)evd.qgpd(x[i], loc=thresh[i], scale=sig[i], shape=xi[i]))
  
  checkEqualsNumeric(eq,myq,msg="qgpd: vectorisation")
  
  #*************************************************************
  # 6.6a Test log.p argument
  
  lq <- qgpd(log(x), sig,xi,thresh,log.p=TRUE)
  
  checkEqualsNumeric(myq, lq, msg="qgpd: log.p=TRUE")
  
  #*************************************************************
  # 6.6a Test log.p argument
  
  LTq <- qgpd(1-x, sig,xi,thresh, lower.tail=FALSE)
  
  checkEqualsNumeric(myq, LTq, msg="qgpd: lower.tail=FALSE")
  
}

test.revTransform <- function(){
  set.seed(20111010)
  n <- 5000
  x <- cbind(rexp(n),rexp(n,3))
  
  x.fit <- migpd(x,mqu = 0.5,penalty="none")
  
  y.l.m <- mexTransform(x.fit,method="mixture", margins="laplace")$transformed
  y.l.e <- mexTransform(x.fit,method="empirical",margins="laplace")$transformed
  
  y.g.m <- mexTransform(x.fit,method="mixture", margins="gumbel")$transformed
  y.g.e <- mexTransform(x.fit,method="empirical",margins="gumbel")$transformed
  
  distFun.l <- function(x) ifelse(x<0, exp(x)/2, 1-exp(-x)/2)
  distFun.g <- function(x) exp(-exp(-x))
  
  u.g.m <- distFun.g(y.g.m)
  u.g.e <- distFun.g(y.g.e)
  
  u.l.m <- distFun.l(y.l.m)
  u.l.e <- distFun.l(y.l.e)
  
  x.l.m <- cbind(revTransform(u.l.m[,1],x[,1],x.fit$mqu[1],x.fit$mth[1],exp(x.fit$models[[1]]$coefficients[1]), x.fit$models[[1]]$coefficients[2]),
                 revTransform(u.l.m[,2],x[,2],x.fit$mqu[2],x.fit$mth[2],exp(x.fit$models[[2]]$coefficients[1]), x.fit$models[[2]]$coefficients[2]))
  x.l.e <- cbind(revTransform(u.l.e[,1],x[,1],method="empirical"),
                 revTransform(u.l.e[,2],x[,2],method="empirical"))
  
  x.g.m <- cbind(revTransform(u.g.m[,1],x[,1],x.fit$mqu[1],x.fit$mth[1],exp(x.fit$models[[1]]$coefficients[1]), x.fit$models[[1]]$coefficients[2]),
                 revTransform(u.g.m[,2],x[,2],x.fit$mqu[2],x.fit$mth[2],exp(x.fit$models[[2]]$coefficients[1]), x.fit$models[[2]]$coefficients[2]))
  x.g.e <- cbind(revTransform(u.g.e[,1],x[,1],method="empirical"),
                 revTransform(u.g.e[,2],x[,2],method="empirical"))
  
  checkEqualsNumeric(x.l.e,x,tolerance=0.0001,msg="revTransform: empirical transformation, laplace target")
  checkEqualsNumeric(x.l.m,x,tolerance=0.0001,msg="revTransform: mixture transformation, laplace target")
  checkEqualsNumeric(x.g.e,x,tolerance=0.0001,msg="revTransform: empirical transformation, gumbel target")
  checkEqualsNumeric(x.g.m,x,tolerance=0.0001,msg="revTransform: mixture transformation, gumbel target")
}

test.rgev <- function() {
  ## so, how do we test an RNG...
  num.simple <- 1000
  num.quantile <- 1e6
  
  xi.values  <- c(0, seq(-5, 5, length.out=10))
  test.quantiles <- c(0.25, 0.5, 0.75)
  
  core.sanity.test <- function(xi) {
    seed <- as.integer(runif(1, -1, 1)*(2**30))
    set.seed(seed)
    samples <- rgev(num.simple, 0, 1, xi)
    checkEquals(length(samples), num.simple,
                "rgev: output of correct length")
    if (xi > 0) {
      checkTrue(all(samples >= -1/xi), "rgev: lower bound check")
    } else if (xi < 0) {
      checkTrue(all(samples <= -1/xi), "rgev: upper bound check")
    }
    ## scale and shift property
    sigma <- rexp(1)
    mu    <- runif(1, -5, 5)
    shifted <- mu + sigma * samples
    set.seed(seed)
    checkEqualsNumeric(shifted,
                       rgev(num.simple, mu, sigma, xi),
                       "rgev: scale and shift")
  }
  
  quantile.test <- function(xi) {
    ## here are the sampled quantiles
    quantiles <- quantile(pgev(rgev(num.quantile, 0, 1, xi),
                               0, 1, xi),
                          probs=test.quantiles,
                          names=FALSE)
    ## this is a bit crude, but hey...
    checkEqualsNumeric(test.quantiles, quantiles,
                       tolerance=0.02,
                       "rgev: quantile test")
  }
  lapply(xi.values, core.sanity.test)
  lapply(xi.values, quantile.test)
}


test.rgpd <- function(){
  ## testing an RNG...
  num.simple <- 1000
  num.quantile <- 1e6
  
  xi.values <- c(0, seq(-2, 2), length.out=10)
  test.quantiles <- c(0.25, 0.5, 0.75)
  
  core.sanity.test <- function(xi) {
    seed <- as.integer(runif(1, -1, 1)*(2**30))
    set.seed(seed)
    samples <- rgpd(num.simple, 1, xi)
    checkEquals(length(samples), num.simple,
                "rgpd: output of correct length")
    if (xi < 0) {
      checkTrue(all(samples <= -1/xi), "rgpd: upper bound check")
    }
    checkTrue(all(samples > 0), "rgpd: lower bound check")
    
    sigma <- rexp(1)
    mu    <- runif(1, -5, 5)
    shifted <- mu + sigma * samples
    set.seed(seed)
    checkEqualsNumeric(shifted,
                       rgpd(num.simple, sigma, xi, u=mu),
                       "rgpd: scale and shift")
  }
  
  quantile.test <- function(xi) {
    ## here are the sampled quantiles
    quantiles <- quantile(pgpd(rgpd(num.quantile, 1, xi),
                               1, xi),
                          probs=test.quantiles,
                          names=FALSE)
    ## this is a bit crude, but hey...
    checkEqualsNumeric(test.quantiles, quantiles,
                       tolerance=0.02,
                       "rgpd: quantile test")
  }
  
  lapply(xi.values, core.sanity.test)
  lapply(xi.values, quantile.test)
}
test.exprel <- function(x) {
  ## pull in from the namespace because RUnit doesn't give access to
  ## internal functions.
  exprel <- texmex:::.exprel
  
  ## first check some simple values
  values <- c(-Inf, NA, seq(-5, 5, length.out=10))
  r.values <- expm1(values) / values
  exprel.values <- exprel(values)
  checkEqualsNumeric(exprel.values, r.values,
                     msg="exprel: value tests")
  
  
  ## the naive R code doesn't compute these right so force the
  ## computation.
  values <- c(Inf, 0)
  exprel.values <- exprel(values)
  desired.values <- c(Inf, 1)
  checkEqualsNumeric(exprel.values, desired.values,
                     msg="exprel: special value tests")
  
}


test.log1prel <- function() {
  ## pull in from the namespace
  
  log1prel <- texmex:::.log1prel
  exprel   <- texmex:::.exprel
  
  ## check simple values
  
  x <- runif(10, -1, 5)
  log1prel.values <- log1prel(x)
  
  r.values <- rep.int(1, length(x))
  flag <- x!=0
  r.values[flag] <- log1p(x[flag]) / x[flag]
  checkEqualsNumeric(log1prel.values, r.values,
                     msg="log1prel: value tests")
  
  ## special values
  checkEqualsNumeric(log1prel(c(-1, 0, Inf)), c(Inf, 1, 0),
                     msg="log1prel: special value tests")
  
  ## check with reference to .exprel
  
  x <- runif(10, -5, 5)
  y <- runif(length(x), -5, 5)
  
  xi <- exprel(x * y) * y
  
  checkEqualsNumeric(log1prel(xi * x) * xi,
                     y, "log1prel: exprel inversion")
  
}


test.specfun.safe.product <- function() {
  prod <- texmex:::.specfun.safe.product
  
  ## simple values
  x <- runif(10, -5, 5)
  y <- runif(length(x), -5, 5)
  
  checkEqualsNumeric(prod(x, y), pmax(x*y, -1),
                     "safe product: simple values")
  
  ## complicated values
  x <- c(0, 0, 1, 1)
  y <- c(Inf, -Inf, Inf, -Inf)
  
  res <- c(0, 0, Inf, -1)
  checkEqualsNumeric(prod(x, y), res,
                     "safe product: complicated values")
}

test.thinAndBurn.evmSim <- function(){
  
  # generate data to use for checking
  d <- sample(3:10,1)
  nrow <- 100
  x <- list(chains = apply(matrix(rep(1:d,each=nrow),ncol=d),2, function( o ) o*1:nrow))
  oldClass( x ) <- "evmSim"
  
  # test appropriate errors for misspecification of thin and burn
  
  checkException(thinAndBurn(x,burn=2),silent=TRUE,msg="thinAndBurn.evmSim: errors for misspecification of thin and burn")
  checkException(thinAndBurn(x,thin=1),silent=TRUE,msg="thinAndBurn.evmSim: errors for misspecification of thin and burn")
  
  #  test burn in
  burn <- sample(nrow/2,1)
  burnOnly <- thinAndBurn(x,burn=burn,thin=1)
  checkEqualsNumeric(x$chains[burn+1,], burnOnly$param[1,],msg="thinAndBurn.evmSim: burn in  ")
  
  # test thinning
  thin <- 2
  thinOnly <- thinAndBurn(x,thin=thin,burn=0)
  
  checkEqualsNumeric(seq(thin,nrow,by=thin), thinOnly$param[,1],msg="thinAndBurn.evmSim: thinning  ")
  
  # test thinning and burning simultaneously
  
  thinBurn <- thinAndBurn(x,thin=thin,burn=burn)
  
  checkEqualsNumeric(seq(burn + thin, nrow,by=thin), thinBurn$param[,1],msg="thinAndBurn.evmSim: thinning and burning simultaneously")
  
  # test returned values of thin and burn
  
  checkEqualsNumeric(thin, thinBurn$thin,burn=0,msg="thinAndBurn.evmSim: test returned value of thin")
  checkEqualsNumeric(burn, thinBurn$burn,thin=1,msg="thinAndBurn.evmSim: test returned value of burn")
  
  # test passing thin and burn via object
  
  x$thin <- thin
  x$burn <- burn
  thinBurn1 <- thinAndBurn(x)
  
  checkEqualsNumeric(seq(burn + thin, nrow,by=thin), thinBurn1$param[,1],msg="thinAndBurn.evmSim: test passing thin and burn via object")
  checkEqualsNumeric(dim(thinBurn$param),dim(thinBurn1$param),msg="thinAndBurn.evmSim: test passing thin and burn via object")
  
  # test thinning and burning a previously thinned and burned object
  
  thin2 <- 4
  burn2 <- 4
  thinBurn2 <- thinAndBurn(thinBurn1,thin=thin2,burn=burn2)
  checkEqualsNumeric(dim(thinBurn1$chains)[2], dim(thinBurn2$chains)[2],msg="thinAndBurn.evmSim: thinning and burning a previously thinned and burned object")
  checkEqualsNumeric((nrow - burn2) / thin2, dim(thinBurn2$param)[1],msg="thinAndBurn.evmSim: thinning and burning a previously thinned and burned object")
}




