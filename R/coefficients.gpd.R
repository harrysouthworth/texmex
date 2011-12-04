`coefficients.gpd` <-
function( object, ... ) {
    res <- object$coefficients
    if ( length( res ) == 2 )
    names( res ) <- c( "phi", "xi" )
    res
}

`coef.gpd` <-
function( object, ... ) {
  coefficients.gpd(object)
}

