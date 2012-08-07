`mexDependence` <-
function (x, which, dqu, margins = "laplace", constrain=TRUE, v = 10, maxit=1000000, start=c(.01, .01), marTransform="mixture", nOptim = 1,
          PlotLikDo=FALSE, PlotLikRange=list(a=c(-1,1),b=c(-3,1)), PlotLikTitle=NULL)
{
   theCall <- match.call()
   if (class(x) != "migpd")
       stop("you need to use an object created by migpd")

   x$margins <-  casefold(margins)
   x <- mexTransform(x, margins = casefold(margins),method = marTransform)

   if (margins == "gumbel" & constrain){
     warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
     constrain <- FALSE
   }

   if (missing(which)) {
       cat("Missing 'which'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
       which <- 1
   }
   else if (length(which) > 1)
       stop("which must be of length 1")
   else if (is.character(which))
       which <- match(which, dimnames(x$transformed)[[2]])

   if (missing(dqu)) {
       cat("Assuming same quantile for thesholding as was used to fit corresponding marginal model...\n")
       dqu <- x$mqu[which]
   }
   dth <- quantile(x$transformed[, which], dqu)

   dependent <- (1:(dim(x$data)[[2]]))[-which]
   if (length(dqu) < length(dependent))
       dqu <- rep(dqu, length = length(dependent))

   # Allowable range of 'a' depends on marginal distributions
   aLow <- ifelse(x$margins == "gumbel", 10^(-10),-1 + 10^(-10))

   if (missing(start)){
     start <- c(.01, .01)
   } else if(class(start) == "mex"){
     start <- start$dependence$coefficients[1:2,]
   }

   if( length(start) == 2 ){
     start <- matrix(rep(start,length(dependent)),nrow=2)
   }

   if( length(start) != 2*length(dependent)){
     stop("start should be of type 'mex' or be a vector of length 2, or be a matrix with 2 rows and ncol equal to the number of dependence models to be estimated")
   }

   if( ! missing(PlotLikRange) ){
     PlotLikDo <- TRUE
   }

   qfun <- function(X, yex, wh, aLow, margins, constrain, v, maxit, start){
     Qpos <- function(param, yex, ydep, constrain, v, aLow) {

  	   a <- param[1]
       b <- param[2]

       res <- PosGumb.Laplace.negProfileLogLik(yex, ydep, a, b, constrain, v, aLow) # defined in file mexDependenceLowLevelFunctions
       res$profLik
     } # Close Qpos <- function

     o <- try(optim(par=start, fn=Qpos,
              control=list(maxit=maxit),
              yex = yex[wh], ydep = X[wh], constrain=constrain, v=v, aLow=aLow), silent=TRUE)

     if (class(o) == "try-error"){
        warning("Error in optim call from mexDependence")
        o <- as.list(o)
        o$par <- rep(NA, 6)
        o$value <- NA
     } else if (o$convergence != 0) {
        warning("Non-convergence in mexDependence")
        o <- as.list(o)
        o$par <- rep(NA, 6)

     } else if(nOptim > 1) {

        for( i in 2:nOptim ){
           o <- try(optim(par=o$par, fn=Qpos,
                    control=list(maxit=maxit),
                    yex = yex[wh], ydep = X[wh], constrain=constrain, v=v, aLow=aLow), silent=TRUE)
           if (class(o) == "try-error"){
             warning("Error in optim call from mexDependence")
             o <- as.list(o)
             o$par <- rep(NA, 6)
             o$value <- NA
             break()
           } else if (o$convergence != 0) {
             warning("Non-convergence in mexDependence")
             o <- as.list(o)
             o$par <- rep(NA, 6)
             break()
           }
        }
     }

     if ( PlotLikDo ){# plot profile likelihood for (a,b)
       nGridPlotLik <- 50
       a.grid <- seq(PlotLikRange$a[1],PlotLikRange$a[2],length=nGridPlotLik)
       b.grid <- seq(PlotLikRange$b[1],PlotLikRange$b[2],length=nGridPlotLik)
       NegProfLik <- matrix(0,nrow=nGridPlotLik,ncol=nGridPlotLik)
       for(i in 1:nGridPlotLik){
         for(j in 1:nGridPlotLik){
           NegProfLik[i,j] <- PosGumb.Laplace.negProfileLogLik(yex=yex[wh], ydep=X[wh],
                                  a = a.grid[i],b=b.grid[j], constrain=constrain,v=v,aLow=aLow)$profLik
         }
       }
       NegProfLik[NegProfLik > 10^10] <- NA
       if(sum(!is.na(NegProfLik))){
          filled.contour(a.grid,b.grid,-NegProfLik,main=paste("Profile likelihood",PlotLikTitle),color.palette = terrain.colors,
                         xlab="a",ylab="b",plot.axes={ axis(1); axis(2); points(o$par[1],o$par[2]) })
       }
     }

     if (!is.na(o$par[1])) { # gumbel margins and negative dependence
        if (margins == "gumbel" & o$par[1] <= 10^(-5) & o$par[2] < 0) {
           lo <- c(10^(-10), -Inf, -Inf, 10^(-10), -Inf, 10^(-10))
           Qneg <- function(yex, ydep, param) {
               param <- param[-1]
               b <- param[1]
               cee <- param[2]
               d <- param[3]
               m <- param[4]
               s <- param[5]

               obj <- function(yex, ydep, b, cee, d, m, s) {
                      mu <- cee - d * log(yex) + m * yex^b
                      sig <- s * yex^b
                      log(sig) + 0.5 * ((ydep - mu)/sig)^2
                      }
               res <- sum(obj(yex, ydep, b, cee, d, m, s))
               res
          }
          o <- try(optim(c(0, 0, 0, 0, 0, 1), Qneg, method = "L-BFGS-B", lower=lo,
                   upper=c(1, 1-10^(-10), Inf, 1-10^(-10), Inf, Inf),
                   yex = yex[wh], ydep = X[wh]), silent=TRUE)

          if (class(o) == "try-error" || o$convergence != 0) {
             warning("Non-convergence in mexDependence")
             o <- as.list(o)
             o$par <- rep(NA, 6)
          }
        } else { # end if gumbel margins and neg dependence
          Z <- (X[wh] - yex[wh] * o$par[1]) / (yex[wh]^o$par[2])
          o$par <- c(o$par[1:2], 0, 0, mean(Z),sd(Z))
        }
    }
    c(o$par[1:6], o$value) # Parameters and negative loglik
   } # Close qfun <- function(

   yex <- c(x$transformed[, which])
   wh <- yex > unique(dth)

   res <- sapply(1:length(dependent),
                 function(X,dat,yex,wh,aLow,margins,constrain,v,maxit,start)qfun(dat[,X],yex,wh,aLow,margins,constrain,v,maxit,start[,X]),
                 dat=as.matrix(x$transformed[, dependent]), yex=yex, wh=wh, aLow=aLow, margins=margins, 
                 constrain=constrain, v=v, maxit=maxit, start=start)

   loglik <- -res[7,]
   res <- matrix(res[1:6,], nrow=6)

   dimnames(res)[[1]] <- c(letters[1:4],"m","s")
   dimnames(res)[[2]] <- dimnames(x$transformed)[[2]][dependent]
   gdata <- as.matrix(x$transformed[wh, -which])
   tfun <- function(i, data, yex, a, b, cee, d) {
       data <- data[, i]
       a <- a[i]
       b <- b[i]
       cee <- cee[i]
       d <- d[i]
       if (is.na(a))
           rep(NA, length(data))
       else {
           if (a < 10^(-5) & b < 0)
               a <- cee - d * log(yex)
           else a <- a * yex
           (data - a)/(yex^b)
       }
   }
   z <- try(sapply(1:(dim(gdata)[[2]]), tfun, data = gdata,
       yex = yex[wh], a = res[1, ], b = res[2, ], cee = res[3, ], d = res[4, ]))
   if (class(z) %in% c("Error", "try-error")) {
       z <- matrix(nrow = 0, ncol = dim(x$data)[[2]] - 1)
   }
   else if (is.R()) {
       if (!is.array(z)) {
           z <- matrix(nrow = 0, ncol = dim(x$data)[[2]] - 1)
       }
   }
   dimnames(z) <- list(NULL,dimnames(x$transformed)[[2]][dependent])
   res2 <- list(coefficients = res, Z = z, dth = unique(dth),
               dqu = unique(dqu), which = which, conditioningVariable= colnames(x$data)[which],
	             loglik=loglik, margins=margins, constrain=constrain, v=v)
   oldClass(res2) <- "mexDependence"

   output <- list(margins=x, dependence=res2, call=theCall)
   oldClass(output) <- "mex"
   output
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

  wavesurge.fit <- migpd(wavesurge,mqu=.7)
  dqu <- 0.8
  which <- 1
  wavesurge.mex <- mexDependence(wavesurge.fit,which=which,dqu=dqu)

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
