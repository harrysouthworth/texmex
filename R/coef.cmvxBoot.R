`coef.cmvxBoot` <-
function( x, which="gpd" ){
	if ( casefold( which ) == "gpd" | which == 1 ) which <- 1
	else which <- 2

	d2 <- dim( x$boot[[ 1 ]][[ which ]] )
	
	B <- x$B
	
    if (which == 1){
        sco <- x$simpleMar[3:4, ] # Point estimates from (penalized) mle
    }
	else {
	    sco <- x$simpleDep # Point estimates of dependence structure
	}
	
	x <- x$boot
	wh <- unlist( lapply( x , function( z ) dim( z$Z )[[ 1 ]] == dim( na.omit( z$Z ) )[[ 1 ]] ) )
	x <- x[ wh ]
	
	eff <- sum( wh )
	
	co <- unlist( lapply( x , function( z, wh ) z[[ wh ]], wh=which ) )
	co <- array( co, dim = c( d2[ 1 ] , d2[ 2 ] , length( co ) / prod( d2 ) ) )
	mco <- apply( co, c( 1, 2 ), mean )
	seco <- apply( co, c( 1 , 2 ), function( x ) sqrt( var( x ) ) )

    # Mean of bootstrap estimates can be used to estimate bias.
    # Then need to subtract bias from parameter estimates.
    b <- mco - sco
    mco <- sco - b

	res <- list( bPoint = mco , se = seco )

	res <- lapply( res , function( x, nms ){
							dimnames( x ) <- nms
							x
						 } ,
						 nms = dimnames( x[[ 1 ]][[ which ]] )
				 )

	attr( res , "Samples generated" ) <- B
	attr( res, "Effective samples" ) <- eff

	res
}

