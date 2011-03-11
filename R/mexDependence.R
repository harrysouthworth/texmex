`mexDependence` <-
function (x, which, dth, dqu, constrain=TRUE, v = 10, maxit=10000)
{
   theCall <- match.call()
   if (class(x) != "migpd")
       stop("you need to use an object created by migpd")

   if (x$margins == "gumbel" & (!missing(constrain) | !missing(v))) {
	   stop("With Gumbel margins, you can't constrain")
   }
   
   if (x$margins == "gumbel"){ constrain <- FALSE }
   
   if (missing(which)) {
       cat("Missing 'which'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
       which <- 1
   }
   else if (length(which) > 1)
       stop("which must be of length 1")
   else if (is.character(which))
       which <- match(which, dimnames(x$transformed)[[2]])
   if (missing(dth) & missing(dqu)) {
       cat("Assuming same quantile for thesholding as was used to fit corresponding marginal model...\n")
       dqu <- x$mqu[which]
   }
   else if (missing(dqu))
       dqu <- x$mqu[which]
   if (missing(dth))
       dth <- quantile(x$transformed[, which], dqu)
   dependent <- (1:(dim(x$data)[[2]]))[-which]
   if (length(dqu) < length(dependent))
       dqu <- rep(dqu, length = length(dependent))

   # Allowable range of 'a' depends on marginal distributions
   lo <- if (x$margins == "gumbel"){ c(10^(-10), -(10^10), -(10^10), 10^(-10)) }
		   else { c(-1 + 10^(-10), -(10^10), -(10^10), 10^(-10)) }
      
   qfun <- function(X, yex, wh, lo, margins, constrain, v, maxit) {
       Q <- function(yex, ydep, param, constrain, v) {

		   a <- param[1]
           b <- param[2]
           m <- param[3]
           s <- param[4]

		   if (a >= 1){ a <- 1 - 10^(-10) }
		   else if (!constrain & a <= 10^(-10)){ a <- 10^(-10) }
		   else if (constrain & a <= -1 + 10^(-10)){ a <- -1 + 10^(-10) }
		   if (b >= 1){ b <- 1 - 10^(-10) }

		   obj <- function(yex, ydep, a, b, m, s, constrain, v) {
               mu <- a * yex + m * yex^b
               sig <- s * yex^b
			   
			   res <- log(sig) + 0.5 * ((ydep - mu)/sig)^2

			   if (constrain){
				   v <- v * max(yex)
				   zpos <- range(ydep - yex)
				   z <- range((ydep - yex * a) / (yex^b)) # q0 & q1 
				   zneg <- range(-ydep - yex) # q0 & q1
				   
				   C1e <- a <= min(1, 1 - b*min(z)*v^(b-1), 1 - v^(b-1)*min(z) + min(zpos)/v) &
						  a <= min(1, 1 - b*max(z)*v^(b-1), 1 - v^(b-1)*max(z) + max(zpos)/v)
				   
				   C1o <- a <= 1 & a > 1 - b * min(z) * v^(b-1) & a > 1 -b*max(z)*v^(b-1) &
							  (1 - 1/b)*(b*min(z))^(1/(1-b)) * (1-a)^(-b/(1 - b)) + min(zpos) > 0 &
							  (1 - 1/b)*(b*max(z))^(1/(1-b)) * (1-a)^(-b/(1 - b)) + max(zpos) > 0

				   C2e <- -a <= min(1, 1 + b*v^(b-1)*min(z), 1 + v^(b-1)*min(z) - min(zneg)/v) &
						  -a <= min(1, 1 + b*v^(b-1)*max(z), 1 + v^(b-1)*max(z) - max(zneg)/v)

				   C2o <- -a <= 1 & -a > 1 + b*v^(b-1)*min(z) & -a > 1 + b*v^(b-1)*max(z) &
						  (1-1/b)*(-b*min(z))^(1/(1-b))*(1+a)^(-b/(1-b)) - min(zneg) > 0 &
						  (1-1/b)*(-b*max(z))^(1/(1-b))*(1+a)^(-b/(1-b)) - max(zneg) > 0

					if (any(is.na(c(C1e, C1o, C2e, C2o)))) {
						warning("Strayed into impossible area of parameter space")
						C1e <- C1o <- C2e <- C2o <- FALSE
					}

				   if (!((C1e | C1o) && (C2e | C2o))){ res <- 10^10 }

			   } # Close if (constrain
			   res
           } # Close obj <- function

		   res <- sum(obj(yex, ydep, a, b, m, s, constrain, v))
               if (is.infinite(res)){
                       if (res < 0){ res <- -(10^10) }
                       else res <- 10^8
                       warning("Infinite value of Q in mexDependence")
               }
           res
       } # Close Q <- function

       o <- try(optim(par=c(0.01, 0.01, 0.01, 1), fn=Q, 
					  method = "Nelder-Mead", control=list(maxit=maxit),
           yex = yex[wh], ydep = X[wh], constrain=constrain, v=v), silent=TRUE)

		if (class(o) == "try-error"){
			warning("Error in optim")
			o <- as.list(o)
			o$par <- rep(NA, 4)
		}
       else if (o$convergence != 0) {
           warning("Non-convergence in mexDependence")
           o <- as.list(o)
           o$par <- rep(NA, 4)
       }
	   
       if (!is.na(o$par[1]))
           if (margins == "gumbel" & o$par[1] <= 10^(-10) & o$par[2] < 0) {
               Q <- function(yex, ydep, param) {
                 param <- param[-1]
                 b <- param[1]
                 cee <- param[2]
                 d <- param[3]
                 m <- param[4]
                 s <- param[5]

				 if (b >=1){ b <- 1 - 10^(-10) }
				 if (d <= 10^(-10)){ d <- 10^(-10) }
				 else if (d >= 1- 10^(-10)){ d <- 1 - 10^(-10) }
			
				 obj <- function(yex, ydep, a, b, cee, d, m, s) {
                   mu <- cee - d * log(yex) + m * yex^b
                   sig <- s * yex^b
                   log(sig) + 0.5 * ((ydep - mu)/sig)^2
                 }
                 res <- sum(obj(yex, ydep, a, b, cee, d, m, s))
                 res
               }
               o <- try(optim(c(0, 0, 0, 0, 0, 1), Q, method = "Nelder-Mead",
                 yex = yex[wh], ydep = X[wh]), silent=TRUE)
               if (class(o) == "try-error" || o$convergence != 0) {
                 warning("Non-convergence in mexDependence")
                 o <- as.list(o)
                 o$par <- rep(NA, 4)
               }
           }
           else o$par <- c(o$par[1:2], 0, 0)
       c(o$par[1:4], o$value) # Parameters and negative loglik
   } # Close qfun <- function(

   
   yex <- c(x$transformed[, which])
   wh <- yex > unique(dth)

	res <- apply(as.matrix(x$transformed[, dependent]), 2, qfun, yex = yex, wh = wh,
	   			lo=lo, margins=x$margins, constrain=constrain, v=v, maxit=maxit)

   loglik <- -res[5,]
   res <- matrix(res[1:4,], nrow=4)

   dimnames(res)[[1]] <- letters[1:4]
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
   res <- list(call = theCall, coefficients = res, Z = z, migpd=x, dth = unique(dth),
       dqu = unique(dqu), which = which, conditioningVariable= colnames(x$data)[which],
	   loglik=loglik)
   oldClass(res) <- "mexDependence"
   res
}

test.mexDependence <- function(){
  smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
  wmarmod <- migpd(winter, mqu=.7,  penalty="none")

  mySdepO3 <- mexDependence(smarmod,which=1,dqu=0.7)
  myWdepO3 <- mexDependence(wmarmod,which=1,dqu=0.7)

  mySdepNO2 <- mexDependence(smarmod,which=2,dqu=0.7)
  myWdepNO2 <- mexDependence(wmarmod,which=2,dqu=0.7)

  mySdepNO <- mexDependence(smarmod,which=3,dqu=0.7)
  myWdepNO <- mexDependence(wmarmod,which=3,dqu=0.7)

  mySdepSO2 <- mexDependence(smarmod,which=4,dqu=0.7)
  myWdepSO2 <- mexDependence(wmarmod,which=4,dqu=0.7)

  mySdepPM10 <- mexDependence(smarmod,which=5,dqu=0.7)
  myWdepPM10 <- mexDependence(wmarmod,which=5,dqu=0.7)

  
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
  
  tol <- 0.11
  if(FALSE){
  par(mfrow=c(2,5))
  plot(jhWdepO3,  myWdepO3$coefficients);abline(0,1)
  plot(jhWdepNO2, myWdepNO2$coefficients);abline(0,1)
  plot(jhWdepNO,  myWdepNO$coefficients);abline(0,1)
  plot(jhWdepSO2, myWdepSO2$coefficients);abline(0,1)
  plot(jhWdepPM10,myWdepPM10$coefficients);abline(0,1)

  plot(jhSdepO3,  mySdepO3$coefficients);abline(0,1)
  plot(jhSdepNO2, mySdepNO2$coefficients);abline(0,1)
  plot(jhSdepNO,  mySdepNO$coefficients);abline(0,1)
  plot(jhSdepSO2, mySdepSO2$coefficients);abline(0,1)
  plot(jhSdepPM10,mySdepPM10$coefficients);abline(0,1)
  }
  
  checkEqualsNumeric(jhWdepO3,  myWdepO3$coefficients,  tol=tol,msg="mexDependence: Winter O3")
  checkEqualsNumeric(jhWdepNO2, myWdepNO2$coefficients, tol=tol,msg="mexDependence: Winter NO2")
  checkEqualsNumeric(jhWdepNO,  myWdepNO$coefficients,  tol=tol,msg="mexDependence: Winter NO")
  checkEqualsNumeric(jhWdepSO2, myWdepSO2$coefficients, tol=tol,msg="mexDependence: Winter SO2")
  checkEqualsNumeric(jhWdepPM10,myWdepPM10$coefficients,tol=tol,msg="mexDependence: Winter PM10")

  checkEqualsNumeric(jhSdepO3,  mySdepO3$coefficients,  tol=tol,msg="mexDependence: Summer O3")
  checkEqualsNumeric(jhSdepNO2, mySdepNO2$coefficients, tol=tol,msg="mexDependence: Summer NO2")
  checkEqualsNumeric(jhSdepNO,  mySdepNO$coefficients,  tol=tol,msg="mexDependence: Summer NO")
  checkEqualsNumeric(jhSdepSO2, mySdepSO2$coefficients, tol=tol,msg="mexDependence: Summer SO2")
  checkEqualsNumeric(jhSdepPM10,mySdepPM10$coefficients,tol=tol,msg="mexDependence: Summer PM10")
  
# test functionality with 2-d data

  wavesurge.fit <- migpd(wavesurge,mq=.7)
  dqu <- 0.8
  which <- 1
  wavesurge.mex <- mexDependence(wavesurge.fit,which=which,dqu=dqu)

  checkEqualsNumeric(dim(wavesurge.mex$Z),c(578,1),msg="mexDependence: execution for 2-d data")
  checkEqualsNumeric(wavesurge.mex$dqu, dqu, msg="mexDependence: execution for 2-d data")
  checkEqualsNumeric(wavesurge.mex$which,which,msg="mexDependence: execution for 2-d data")
  
}
