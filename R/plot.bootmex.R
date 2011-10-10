`plot.bootmex` <-
function( x , plots = "gpd", main="", ... ){

    # Want to look at the marginal GPD parameters or the
    # dependence structure parameters?
	if ( casefold( plots ) == "gpd" ) { which <- 1 }
	else { which <- 2 }
	
	d2 <- dim(x$boot[[1]][[which]])
	pointEst <- x$simpleDep
  
  condVar <- names(x$simpleMar$data)[x$which]
  margins <- x$margins
	x <- x$boot
	co <- unlist( lapply( x , function( z, wh ) z[[ wh ]], wh=which ) )
	co <- array(co, dim = c( d2[1] , d2[2] , length(co) / prod(d2)))

  lco <- list(length=prod(d2))

	for (i in 1:d2[2]){ # loop over variables
	  for (j in 1:d2[1]){ # loop over parameters
	    lco[[ j + d2[1]*(i - 1) ]] <- co[j, i, ]
	  } # close j
  } # close i
	
  cn <- colnames(x[[1]][[which]]) # variable names
  rn <- rownames(x[[1]][[which]]) # parameter names
  if(which == 2){
    cn <- paste(cn, "|", condVar)
  }
  labs <- paste(rep(rn, length(cn)), rep(cn, each=switch(which,2,6)),sep="  ")

	fun <- function(X, z, label, ...) {
		hist(z[[X]] , prob=TRUE, xlab=label[X], main=main, ...)
		lines(density( z[[X]], n=100 ))
		invisible()
	}

  if(which == 1){
    lapply(1:prod(d2), fun, z=lco, label=labs, ...)
  }

  if(which == 2){ # scatterplots of dependence parameters    
    fun <- function(X,z,label, ...){
      offset <- (X-1) * 6
      plot(lco[[offset + 1]],lco[[offset + 2]],xlab=labs[offset + 1],ylab=labs[offset + 2],main=main, ...)
      points(pointEst[1,X],pointEst[2,X],pch="@",col="red")
      if( margins == "gumbel"){
        plot(lco[[offset + 3]],lco[[offset + 4]],xlab=labs[offset + 3],ylab=labs[offset + 4],main=main, ...)
        points(pointEst[3,X],pointEst[4,X],pch="@",col="red")
      }
    }
    lapply(1:d2[2], fun, z=lco,label=labs, ...)
  }
	invisible()
}

test.plot.bootmex <- function(){

set.seed(3141593)

# 2-d wavesurge data

  wavesurge.fit <- mex(wavesurge,which=1,mqu=0.7) 
  wavesurge.boot <- bootmex(wavesurge.fit,R=50)
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

  Sboot.O3 <- bootmex(smarmex.O3, R=R)
  Wboot.O3 <- bootmex(wmarmex.O3, R=R)
  Sboot.NO2 <- bootmex(smarmex.NO2, R=R)
  Wboot.NO2 <- bootmex(wmarmex.NO2, R=R)
  Sboot.NO <- bootmex(smarmex.NO, R=R)
  Wboot.NO <- bootmex(wmarmex.NO, R=R)
  Sboot.SO2 <- bootmex(smarmex.SO2, R=R)
  Wboot.SO2 <- bootmex(wmarmex.SO2, R=R)
  Sboot.PM10 <- bootmex(smarmex.PM10, R=R)
  Wboot.PM10 <- bootmex(wmarmex.PM10, R=R)

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
