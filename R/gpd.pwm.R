`gpd.pwm` <-
function( x ){
	m = mean( x )
	x = sort( x )
	n = length( x )
	p = ( 1:n - .35 ) / n
	s = mean( ( 1 - p ) * x )
	xi = ( m / (m - 2 * s ) ) -2
	b = 2 * m * s / ( m - 2 * s )
	if ( x[ n ] < b / xi ){
		res = c( sigma = b, xi = xi )
		attr( res , "type" ) = "PWM"
	}
	else {
		res = c( sigma = b, xi = b/x[ n ] )
		attr( res, "type" ) = "Hybrid"
	}
	res
}

