`show.summary.cmvxPrediction` <-
function( x ){
	print( x$call )
		cat( "\nResults from", length( x$B ), "bootstrap runs.\n" )
	cat( paste( "\nConditioned on ", x$which, " being above its ", 100*x$pqu, "th percentile.\n\n", sep = "" ) )
	cat( "\nExpected values:\n\n" )
	print( x$ans )

	cat( "\nConditional probability of threshold exceedance:\n\n" )

	print( x$thres )
	
	invisible()
}

