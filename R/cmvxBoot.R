`cmvxBoot` <-
function( x, which , B=100, gth, gqu, nPass = 3, trace=10 ){
	theCall <- match.call()

	if ( class( x ) != "migpd" ) stop( "object should have class migpd" )

   if (missing(which)) {
       cat("Missing 'which'. Conditioning on", dimnames(x$gumbel)[[2]][1],
           "\n")
       which <- 1
   }

    if (missing(gqu) & missing(gth)){
        gqu <- x$qu[which]
    }
    if ( missing( gth ) ) gth <- quantile( x$gumbel[, which ] , prob=gqu )


	penalty <- x$penalty
	priorParameters <- x$priorParameters

	# This is brutal and could be made more efficient
	innerFun <- function( i, x, which, gth, gqu, penalty, priorParameters, pass=1, trace=trace ){

		# This could be improved by pushing some of the calculations outside the 
		# function, and by allowing starting values to be passed to the functions
		# that use optim

		n <- dim( x$gumbel )[[ 1 ]]
		d <- dim( x$gumbel )[[ 2 ]]
		
		if ( is.character( which ) )
			which <- ( 1:d )[ dimnames( x$data )[[ 2 ]] == which ]
		dependent <- ( 1:d )[ -which ]
		
		# Take bs sample from Gumbel data and get ranks
		g <- sample( 1:(dim( x$gumbel )[[ 1 ]] ), size=n, replace = TRUE )
		g <- x$gumbel[ g, ]

		ok <- FALSE
        while( !ok ){ # Occasionally get a sample with no extremes
            for( j in 1:( dim( g )[[ 2 ]] ) ) # -log( -log(u) ) ~ G( 0, 1 )
                g[ order( g[,j] ) , j ] <- sort( -log( -log( runif( dim( g )[[ 1 ]] ) ) ) )
            if ( sum( g[ , which ] > gth ) > 1 ) ok <- TRUE
        } # Close while( !ok ) 

		
		# Need to transform to original scale
		getgum <- function( i , x, data, mod, th, qu ){
			x <- c( x[, i ] )
			param <- mod[[ i ]]$par
			th <- th[ i ]
			qu <- qu[ i ]
			data <- c( data[, i ] )

			res <- revGumbel( exp( -exp( -x ) ) , data=data,
								th=th, qu=qu,
								sigma=exp( param[ 1 ] ), xi = param[ 2 ] )
			res
		}

		g <- sapply( 1:d , getgum, x=g, data=x$data, mod=x$models, th=x$th, qu=x$qu )
		dimnames( g )[[ 2 ]] <- names( x$models )
	
		# Get GPD parameters
		ggpd <- migpd( g , th = x$th, penalty=penalty, priorParameters=priorParameters ) # 

		if ( !missing( gqu ) ) gqu <- rep( gqu , length=d )
		# Get dependence parameters
		gth <- unlist( lapply( 1:length( gqu ),
							   function( i, d, qu )
							   	quantile( d[,i ], qu[i] ) ,
							   d=x$gumbel, qu=gqu ) )
		gth <- quantile( c( x$gumbel[, which ] ), gqu[ which ] )
		# Need to pass gth, not gqu, because gth is fixed
		gd <- cmvxDependence( ggpd , gth=gth, which=which )

		# Need to return GPD parameters and the dependence structure parameters
		res <- list( GPD=coef( ggpd )[ 3:4,] , dependence=gd$parameters, Z = gd$Z )

		if ( pass == 1 ){
			if ( i %% trace == 0 ){
				cat( paste( i, "replicates done\n" ) )
                        }
                } # Close if ( pass == 1
		
		res
	} # Close innerFun
	res <- lapply( 1:B, innerFun, x=x, which=which, gth=gth, gqu=gqu,
				   penalty=penalty, priorParameters=priorParameters,
				pass = 1, trace=trace )

	# Re run for non-converged sets
	if ( nPass > 1 ){
		for( pass in 2:nPass ){
			rerun <- sapply( res , function( x ) any( sapply( x , function( x ) any( is.na( x ) ) ) ) )
			wh <- !unlist( lapply( res , function( x ) dim( x$Z )[[ 1 ]] > 0 ) )
			rerun <- apply( cbind( rerun, wh ) , 1, any )
			if ( sum( rerun ) > 0 ){
				cat( "Pass", pass, ":", sum( rerun ), "samples to rerun.\n" )
				rerun <- ( 1:B )[ rerun ]
				res[ rerun ] <- lapply( (1:B)[ rerun ], innerFun, x=x, which=which, gth=gth, gqu=gqu,
						   penalty=penalty, priorParameters=priorParameters , pass=pass)
			}
		} # close for( pass in...
	} # close if (nPass > 1 )
	
#	wh <- unlist( lapply( res , function( x ) dim( x$Z )[[ 1 ]] > 0 ) )
#	res <- res[ wh ]

	ans <- list()
	ans$boot <- res
	ans$call <- theCall
	ans$gqu <- gqu
	ans$which <- which
	ans$B <- B
	ans$simpleEsts <- coef(x)
	ans$simpleDep <- cmvxDependence(x, which)$parameters
	oldClass( ans ) <- "cmvxBoot"
	ans
}

