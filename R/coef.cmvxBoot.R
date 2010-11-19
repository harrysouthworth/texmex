`coef.cmvxBoot` <-
function( x ){
	d2 <- dim( x$boot[[ 1 ]][[ 2 ]] )
	
	B <- x$B
	
  sco <- x$simpleDep # Point estimates of dependence structure
	
	x <- x$boot
	wh <- unlist( lapply( x , function( z ) dim( z$Z )[[ 1 ]] == dim( na.omit( z$Z ) )[[ 1 ]] ) )
	x <- x[ wh ]
	
	eff <- sum( wh )
	
	co <- unlist( lapply( x , function( z ) z[[ 2 ]]) )
	co <- array( co, dim = c( d2[ 1 ] , d2[ 2 ] , length( co ) / prod( d2 ) ) )
	mco <- apply( co, c( 1, 2 ), mean )
	seco <- apply( co, c( 1 , 2 ), function( x ) sqrt( var( x ) ) )

	res <- list( bPoint = mco , se = seco )

	res <- lapply( res , function( x, nms ){
							dimnames( x ) <- nms
							x
						 } ,
						 nms = dimnames( x[[ 1 ]][[ 2 ]] )
				 )

	attr( res , "Samples generated" ) <- B
	attr( res, "Effective samples" ) <- eff

	res
}

