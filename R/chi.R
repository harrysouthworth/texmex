chi <- 
  # Much of the code in here was written by Alec Stephenson
  # and comes from his 'chiplot' function in his 'evd' package.
  # Minor differences between the evd implementation and this are: 
  # use of ties.method="first" here as oppsed to the evd method 
  # which used ties.method="average"; lower bound on chibar wrong 
  # in evd package.
function (data, nq = 100, qlim = NULL, alpha = 0.05, trunc = TRUE) {
    
    theCall <- match.call()

    eps <- .Machine$double.eps^0.5

    data <- na.omit(data)
    n <- nrow(data)

    # Get the EDFs
    t.method <- "first"

	if (is.R()){
	    data <- cbind(rank(data[, 1],ties.method = t.method)/(n + 1), 
    	              rank(data[, 2],ties.method = t.method)/(n + 1))
	}
	else {
		data <- cbind(rank(data[, 1])/(n + 1), 
    	              rank(data[, 2])/(n + 1))
	}

    rowmax <- apply(data, 1, max)
    rowmin <- apply(data, 1, min)

    qlim2 <- c(min(rowmax) + eps, max(rowmin) - eps)

    if (!is.null(qlim)) {
        if (qlim[1] < qlim2[1]){ stop("lower quantile limit is too low") }
        if (qlim[2] > qlim2[2]){ stop("upper quantile limit is too high") }
        if (qlim[1] > qlim[2]){ stop("lower quantile limit is less than upper quantile limit") }
    }
    else{
        qlim <- qlim2
    }

    u <- seq(qlim[1], qlim[2], length = nq)

    # HS. Replaced 2 for loops with calls to sapply

    cu <- sapply(1:nq, function(i, x, y){ mean(y < x[i]) }, y=rowmax, x=u )
    cbaru <- sapply(1:nq, function(i, x, y){ mean(y > x[i]) }, y=rowmin, x=u )

    # Get \chi and \bar\chi
    chiu <- 2 - log(cu)/log(u) # 3.2 of Coles, Heffernan, Tawn
    chibaru <- (2 * log(1 - u))/log(cbaru) - 1 # Page 348 of Coles, Heffernan, Tawn

    # Get confidence limits
    varchi <- ((1/log(u)^2 * 1)/cu^2 * cu * (1 - cu))/n
    varchi <- qnorm(1 - alpha/2) * sqrt(varchi)

    varchibar <- (((4 * log(1 - u)^2)/(log(cbaru)^4 * cbaru^2)) * cbaru * (1 - cbaru))/n
    varchibar <- qnorm(1 - alpha/2) * sqrt(varchibar)

    chiu <- cbind(chilow = chiu - varchi,           # Lower
                  chi = chiu,                       # Point est.
                  chiupp = chiu + varchi)           # Upper

    chibaru <- cbind(chiblow = chibaru - varchibar, # Lower
                     chib = chibaru,                # Point est.
                     chibupp = chibaru + varchibar) # Upper

    chiulb <- 2 - log(pmax(2 * u - 1, 0))/log(u)

    if (trunc) {
        chiu[chiu > 1] <- 1
        chiu <- apply(chiu, 2, function(x, z){ pmax(x, z) }, z = chiulb)
        chibaru[chibaru > 1] <- 1
        chibaru[chibaru < -1] <- -1
    }

    res <- list(chi=chiu, chibar = chibaru, quantile=u, call=theCall, 
                qlim=qlim, chiulb = chiulb)
    oldClass(res) <- "chi"
    res
} # Close chi <- function


print.chi <- function(x, ...){
    print(x$call,...)
    cat("Values of chi and chi-bar obtained and",
         length(x$quantile), "quantiles.\n")
    invisible()
}

summary.chi <- function(object, digits=3, ...){
    print(object$call)
    cat("Values of chi and chi-bar obtained and",
         length(object$quantile), "quantiles.\n")

    wh <- quantile(object$quantile, prob=c(.05, .5, .95))
    wh <- sapply(wh, function(i, u){
                        d <- abs(u - i)
                        min(u[d == min(d)])
                     }, u=object$quantile)

    chiQ <- object$chi[object$quantile %in% wh, 2]
    chibarQ <- object$chibar[object$quantile %in% wh, 2]

    res <- rbind(wh, chiQ, chibarQ)
    dimnames(res) <- list(c("quantile", "chi", "chi-bar"), rep("", 3))
    print(res, digits=digits, ...)
    invisible(res)
}


plot.chi <- function(x, which=1:2, lty = 1, cilty = 2, col = 1, spcases = FALSE, cicol = 1,
                     xlim = c(0, 1), ylim1 = c(-1, 1), ylim2 = c(-1, 1),
                     main1 = "Chi", main2 = "Chi Bar",
                     xlab = "Quantile", 
                     ylab1 = expression(chi),#"Chi",
                     ylab2 = expression(bar(chi)), #"Chi Bar",
                     ask, ...){

    show <- logical(2)
    show[which] <- TRUE
    lty <- c(cilty, lty, cilty)
    col <- c(cicol, col, cicol)
    nb.fig <- prod(par("mfcol"))

	if (is.R() & missing(ask)){
		ask <- nb.fig < length(which) && dev.interactive()
	}
	else {
		ask <- FALSE
	}
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
        matplot(x$quantile, x$chi, type = "l", lty = lty, col = col, xlim = xlim, 
            ylim = ylim1, main = main1, xlab = xlab, ylab = ylab1, 
            ...)
        if (spcases) {
            segments(x$qlim[1], 0, x$qlim[2], 0, lty = 5, col = "grey")
            segments(x$qlim[1], 1, x$qlim[2], 1, lty = 5, col = "grey")
            lines(x$quantile, x$chiulb, lty = 5, col = "grey")
        }
    }
    if (show[2]) {
        matplot(x$quantile, x$chibar, type = "l", lty = lty, col = col, 
            xlim = xlim, ylim = ylim2, main = main2, xlab = xlab, 
            ylab = ylab2, ...)
        if (spcases) {
            segments(x$qlim[1], 0, x$qlim[2], 0, lty = 5, col = "grey")
            segments(x$qlim[1], 1, x$qlim[2], 1, lty = 5, col = "grey")
            lines(x$quantile, x$chibarulb, lty = 5, col = "grey")
        }
    }
    invisible()
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
  res <- plot(chi,main1="Figure 8.11 of Coles (2001)\nChi")
  checkEquals(res,NULL,msg = "plot.chi: check successful execution")
}
