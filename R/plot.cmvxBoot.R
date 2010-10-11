`plot.cmvxBoot` <-
function( x , which = "gpd", ... ){

    # Want to look at the marginal GPD parameters or the
    # dependence structure parameters?
	if ( casefold( which ) == "gpd" ) { which <- 1 }
	else { which <- 2 }
	
	d2 <- dim(x$boot[[1]][[which]])
	
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
    labs <- paste(rep(cn, each=2), rep(rn, length(cn)))




	fun <- function(X, z, label, ...) {
		hist(z[[X]] , prob=TRUE, xlab=label[X], main="", ...)
		lines(density( z[[X]], n=100 ))
		invisible()
	}

			
    wh <- lapply(1:prod(d2), fun, z=lco, label=labs, ...)
#	wh <- apply( co , 1:2, function(x) x )

	invisible()
}
