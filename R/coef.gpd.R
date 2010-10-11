`coef.gpd` <-
function( x ) {
								res = x$mle
								if ( length( res ) == 2 )
									names( res ) = c( "log(sigma)", "xi" )
								res
							   }

