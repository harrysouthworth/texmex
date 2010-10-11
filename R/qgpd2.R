`qgpd2` <-
function(p , sigma = 1 , xi = 1 , u = 0, la = 1 )
	u + ( sigma * ( ( p * la )^( xi ) - 1)) / xi

