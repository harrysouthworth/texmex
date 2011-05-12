`coefficients.gpd` <-
function( object, ... ) {
								res <- object$coefficients
								if ( length( res ) == 2 )
									names( res ) <- c( "log(sigma)", "xi" )
								res
							   }

`coef.gpd` <-
function( object, ... ) {
  coefficients.gpd(object)
							   }

