`coef.mexBoot` <-
function( object, ... ){
	d2 <- dim( object$boot[[ 1 ]][[ 2 ]] )
	
	B <- object$B
	
  sco <- object$simpleDep # Point estimates of dependence structure
	
	object <- object$boot
	wh <- unlist( lapply( object , function( z ) dim( z$Z )[[ 1 ]] == dim( na.omit( z$Z ) )[[ 1 ]] ) )
	object <- object[ wh ]
	
	eff <- sum( wh )
	
	co <- unlist( lapply( object , function( z ) z[[ 2 ]]) )
	co <- array( co, dim = c( d2[ 1 ] , d2[ 2 ] , length( co ) / prod( d2 ) ) )
	mco <- apply( co, c( 1, 2 ), mean )
	seco <- apply( co, c( 1 , 2 ), function( object ) sqrt( var( object ) ) )

	res <- list( bPoint = mco , se = seco )

	res <- lapply( res , function( object, nms ){
							dimnames( object ) <- nms
							object
						 } ,
						 nms = dimnames( object[[ 1 ]][[ 2 ]] )
				 )

	attr( res , "Samples generated" ) <- B
	attr( res, "Effective samples" ) <- eff

	res
}

