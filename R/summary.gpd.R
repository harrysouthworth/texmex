`summary.gpd` <-
function( object , nsim = 1000 , alpha = .050, ... ){
	env <- qqgpd( object, plot = FALSE, nsim = nsim, alpha = alpha )
	res <- list( model = object, envelope = env , nsim = nsim, alpha = alpha )
	oldClass( res ) <- "summary.gpd"
	res
}

