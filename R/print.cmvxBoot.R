`print.cmvxBoot` <-
function( x , ... ){
	print( x$call )
	cat( paste( "\n", x$B, " bootstrap samples created.\n\n" , sep = "" ) )

  co <- coef( x )
  eff <- attributes( co )$"Effective samples"
  cat( paste("Dependence structure estimation successful for", eff, "effective samples.\n" ) )
  cat( paste("Dependence structure bootstrap mean parameter estimates:\n" ))
  print( co[[ 1 ]] )
  cat("\n")

	}

