`plot.bootmex` <-
function( x , plots = "gpd", main="", ... ){

    # Want to look at the marginal GPD parameters or the
    # dependence structure parameters?
	if ( casefold( plots ) == "gpd" ) { which <- 1 }
	else { which <- 2 }
	
	d2 <- dim(x$boot[[1]][[which]])
	pointEst <- x$simpleDep
  
  condVar <- names(x$simpleMar$data)[x$which]
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
  labs <- paste(rep(rn, length(cn)), rep(cn, each=which*2),sep="  ")

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
      offset <- (X-1) * 4
      plot(lco[[offset + 1]],lco[[offset + 2]],xlab=labs[offset + 1],ylab=labs[offset + 2],main=main, ...)
      points(pointEst[1,X],pointEst[2,X],pch="@",col="red")
      plot(lco[[offset + 3]],lco[[offset + 4]],xlab=labs[offset + 3],ylab=labs[offset + 4],main=main, ...)
      points(pointEst[3,X],pointEst[4,X],pch="@",col="red")
    }
    lapply(1:d2[2], fun, z=lco,label=labs, ...)
  }
	invisible()
}

test.plot.bootmex <- function(){

# 2-d wavesurge data

  wavesurge.fit <- migpd(wavesurge,mqu=0.7) 
  wavesurge.boot <- bootmex(wavesurge.fit,which=1,R=50)
  par(mfrow=c(3,2),pty="m")
  check1 <- plot(wavesurge.boot,main="Marginal parameters\nWave surge data of Coles 2001")
  check2 <- plot(wavesurge.boot,plots="dep",main="Dependence parameters\nWave surge data of Coles 2001")
  
# 5-d air pollution data

  smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
  wmarmod <- migpd(winter, mqu=.7,  penalty="none")
  
  Qu <- 0.7
  R <- 50
  Sboot.O3 <- bootmex(smarmod, which=1, dqu=Qu, R=R)
  Wboot.O3 <- bootmex(wmarmod, which=1, dqu=Qu, R=R)
  Sboot.NO2 <- bootmex(smarmod, which=2, dqu=Qu, R=R)
  Wboot.NO2 <- bootmex(wmarmod, which=2, dqu=Qu, R=R)
  Sboot.NO <- bootmex(smarmod, which=3, dqu=Qu, R=R)
  Wboot.NO <- bootmex(wmarmod, which=3, dqu=Qu, R=R)
  Sboot.SO2 <- bootmex(smarmod, which=4, dqu=Qu, R=R)
  Wboot.SO2 <- bootmex(wmarmod, which=4, dqu=Qu, R=R)
  Sboot.PM10 <- bootmex(smarmod, which=5, dqu=Qu, R=R)
  Wboot.PM10 <- bootmex(wmarmod, which=5, dqu=Qu, R=R)

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
