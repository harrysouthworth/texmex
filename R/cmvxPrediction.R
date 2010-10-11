`cmvxPrediction` <-
function( migpd , boot, pqu = .99, nsim = 1000 ){
	theCall <- match.call()
	
	if ( class( migpd ) != "migpd" )
		stop( "you need to use an object with class 'migpd'" )
	if ( class( boot ) != "cmvxBoot" )
		stop( "you need to provide a cmvxBoot object" )
		
	which <- boot$which
	if ( is.character( which ) )
		which <- match( which, dimnames( migpd$data )[[ 2 ]] )
	
	################################################################
	# The function lfun does most of the work
	lfun <- function( i , bo, pqu, nsim , migpd, which ){
		if ( i %% 10 == 0 ) cat( i, "sets done\n" )

		# Step 1 - get y 
		y <- -log( -log( runif( nsim , min=pqu ) ) )
		
		# Step 2 - get z
		z <- bo[[ i ]]$Z
		z <- z[ sample( 1:( dim( z )[[ 1 ]] ), size = nsim, replace = TRUE ), ]

		# Step 3 - get y_{-i}
		ab <- bo[[ i ]]$dependence

		fun <- function( i, z, v , y ){
			v <- v[ , i ]
			z <- z[ , i ]
			if ( !is.na( v[ 1 ] ) ){
				if( v[ 1 ] < 10^(-5) & v[ 2 ] < 0 ){
					if( v[ 4 ] < 10^(-5 ) ) d <- 0
					else d <- v[ 4 ]
					a <- v[ 3 ] - d * log( y )
				}
				else a <- v[ 1 ] * y
			} # close if( !is.na...
#			else a <- NA
			a + ( y^v[ 2 ] ) * z
		}
		ymi <- sapply( 1:( dim( z )[[ 2 ]] ) , fun, z=z, v=ab , y=y )

		# Step 4 - transform to original scale
		# Transform to (0, 1)
		# Do xi
		xp <- pqu * ( 1 - migpd$qu[ which ] )
		cox <- coef( migpd )[, which ]
		m <- 1 / ( 1 - pqu ) # Need to estimate qpu quantile
		zeta <- 1 - migpd$qu[ which ] # Coles, page 81
		xth <- migpd$th[ which ] + cox[ 3 ] / cox[ 4 ] * ( ( m*zeta )^cox[ 4 ] - 1 )
		xi <- u2gpd( exp( -exp( -y ) ), p = 1 - pqu, th=xth, sigma=cox[ 3 ], xi = cox[ 4 ] )

		xmi <- apply( ymi , 2, function( x ) exp( -exp( -x ) ) )

		# Push through CDF to get to original scale
		# We're above the original threshold, so use gpd for x
		for( j in 1:( dim( z )[[ 2 ]] ) ){
			cox <- bo[[ i ]]$GPD[,-which][, j ]
			xmi[, j ] <- revGumbel( xmi[, j ], migpd$data[,-which][, j ],
								  th=migpd$th[ -which ][ j ],
								  qu = migpd$qu[-which][ j ],
								  sigma=cox[ 1 ], xi=cox[ 2 ] )
		}# Close for j
		dimnames( xmi )[[ 2 ]] <- dimnames( coef( migpd ) )[[ 2 ]][ -which ]
		res <- data.frame( xi , xmi )
		names( res ) <- c( dimnames( coef( migpd ) )[[ 2 ]][ which ] , dimnames( xmi )[[ 2 ]] )
		res
	} # Close lfun

	bootRes <- lapply( 1:length( boot$boot ) , lfun ,
				   migpd=migpd, pqu=pqu, bo = boot$boot, nsim=nsim,
				   which = which )
	
	# bootRes contains the bootstrap simulated Y_{-i}

	##########################################################################
	# Also get a sample using the point estimates of the parameters
	# that are suggested by the data
	y <- -log( -log( runif( nsim , min=pqu ) ) )
		
	# Step 2 - get z
	dall <- cmvxDependence( migpd , which=which , gqu=boot$gqu )
	dco <- dall$parameters
	z <- dall$Z
	# Resample rows of Z
	z <- z[ sample( 1:( dim( z )[ 1 ] ), size=nsim, replace=TRUE ) ,]

	# Step 3 - get y_{-i}
	fun <- function( i, z, v , y ){
		v <- v[ , i ]
		z <- z[ , i ]
		if( v[ 1 ] < 10^(-5) & v[ 2 ] < 0 ){
			if( v[ 4 ] < 10^(-5 ) ) d <- 0
			else d <- v[ 4 ]
			a <- v[ 3 ] - d * log( y )
		}
		else a <- v[ 1 ] * y
		a + ( y^v[ 2 ] ) * z
	}
	ymi <- sapply( 1:( dim( z )[[ 2 ]] ) , fun, z=z, v=dco , y=y )
	# Checked dim( ymi ) - ok

	# Step 4 - transform to original scale
	# Transform to (0, 1)
	xi <- exp( -exp( -y ) )
	xmi <- apply( ymi, 2, function( x ) exp( -exp( -x ) ) )

	# Do xi
	xp <- pqu * ( 1 - migpd$qu[ which ] )
	cox <- coef( migpd )[, which ]
	m <- 1 / ( 1 - pqu ) # Need to estimate qpu quantile
	zeta <- 1 - migpd$qu[ which ] # Coles, page 81
	xth <- migpd$th[ which ] + cox[ 3 ] / cox[ 4 ] * ( ( m*zeta )^cox[ 4 ] - 1 )
	xi <- u2gpd( xi, p = 1 - pqu, th=xth, sigma=cox[ 3 ], xi = cox[ 4 ] )
	
	# Now need to transform each column of xmi
	for( i in 1:( dim( xmi )[[ 2 ]] ) ){
		cox <- coef( migpd )[ , -which ][ , i ]
		xmi[, i ] <- revGumbel( xmi[ ,i ], migpd$data[,-which][, i ],
								  th=migpd$th[ -which ][ i ],
								  qu = migpd$qu[-which][ i ],
								  sigma=cox[ 3 ], xi=cox[ 4 ] )
	}

	sim <- data.frame( xi , xmi )
	names( sim ) <- c( names( migpd$data )[ which ],
					   names( migpd$data )[ -which ]
					  )
	
	# Wrap up and return to the user
	data <- list( real = data.frame( migpd$data[, which ], migpd$data[, -which] ) ,
				  simulated = sim, pth=xth
				 )

	res <- list( call = theCall , replicates = bootRes, data = data,
				 which = which, pqu = pqu,
				 th=c( migpd$th[ which ], migpd$th[ -which ] ) )
	
	oldClass( res ) <- "cmvxPrediction"

	invisible( res )
}

