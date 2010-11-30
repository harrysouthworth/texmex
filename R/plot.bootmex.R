`plot.bootmex` <-
function( x , which = "gpd", main="", ... ){

    # Want to look at the marginal GPD parameters or the
    # dependence structure parameters?
	if ( casefold( which ) == "gpd" ) { which <- 1 }
	else { which <- 2 }
	
	d2 <- dim(x$boot[[1]][[which]])
	pointEst <- x$simpleDep
  
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
    labs <- paste(rep(cn, each=4), rep(rn, length(cn)))

	fun <- function(X, z, label, ...) {
		hist(z[[X]] , prob=TRUE, xlab=label[X], main=main, ...)
		lines(density( z[[X]], n=100 ))
		invisible()
	}

  wh <- lapply(1:prod(d2), fun, z=lco, label=labs, ...)

  if(which == 2){ # scatterplots of dependence parameters    
    fun <- function(X,z,label, ...){
      offset <- (X-1) * 4
      plot(lco[[offset + 1]],lco[[offset + 2]],xlab=labs[offset + 1],ylab=labs[offset + 2],main=main, ...)
      points(pointEst[1,X],pointEst[2,X],pch="@",col="red")
      plot(lco[[offset + 3]],lco[[offset + 4]],xlab=labs[offset + 3],ylab=labs[offset + 4],main=main, ...)
      points(pointEst[3,X],pointEst[4,X],pch="@",col="red")
    }
    lapply(1:d2[2], fun, z=lco,label=labs)
  }
	invisible()
}
