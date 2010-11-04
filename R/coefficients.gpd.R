`coefficients.gpd` <-
function( x ) {
								res = x$coefficients
								if ( length( res ) == 2 )
									names( res ) = c( "log(sigma)", "xi" )
								res
							   }

