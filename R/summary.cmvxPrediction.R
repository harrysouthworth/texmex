`summary.cmvxPrediction` <-
function( object, th, probs=c( .05, .5, .95 ), ... ){

	if ( is.R() ) stdev <- function( x ) sqrt( var( x ) )
	if ( missing( th ) ) th <- object$th

	res <- t( sapply( object$replicates , function ( x ) apply( x, 2, mean ) ) )

	sumfun <- function( x , probs){
		c( mean=mean( x ), se=stdev( x ) , quantile( x, probs=probs ) )
    }

	ans <- apply( res, 2, sumfun, probs ) # Summary of expected values
	dn <- paste( "E(", dimnames( object$replicates[[ 1 ]] )[[ 2 ]] ,"|", names( object$data$simulated )[[ 1 ]] , ">Q",100*object$pqu,")", sep="" )
	dimnames( ans )[[ 2 ]] <- dn

	thres <- t( sapply( 1:( dim( object$replicates[[ 1 ]] )[[ 2 ]] ) ,
						function( i, x, th ){
							x <- sapply( x , function( x, i ) x[,i], i=i )
							th <- th[ i ]
							apply( x , 2, function( x, th )
											mean( x > th ),
										th = th )
						}, x=object$replicates, th = th ) )
	thres <- apply( thres, 1, mean )
	thres <- matrix( thres, nrow=1 )

	wn <- dimnames( object$data$simulated )[[ 2 ]][ 1 ]
	wth <- paste( "Q", 100*object$pqu, sep = "" )
#	dn <- paste( "P(", dimnames( object$replicates[[ 1 ]] )[[ 2 ]] , ">", th[ -1 ],"|", wn, ">", wth, ")", sep = "" )
	dn <- paste( "P(", dimnames( object$replicates[[ 1 ]] )[[ 2 ]] , ">", th,"|", wn, ">", wth, ")", sep = "" )

	dimnames( thres ) <- list( "", dn )
	
	ans <- list( ans=ans, thres=thres, call=object$call, pqu=object$pqu ,
				 B = length( object$replicates ),
				 which = names( object$data$simulated )[[ 1 ]],
				 statistic=deparse( substitute( statistic ) )
				 )

	oldClass( ans ) <- "summary.cmvxPrediction"
	ans
}

