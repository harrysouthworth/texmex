`print.bootmex` <-
function( y , ... ){
	print( y$call )
	cat( paste( "\n", y$R, " bootstrap samples created.\n\n" , sep = "" ) )

  eff <- attributes( co )$"Effective samples"
  cat( paste("Dependence structure estimation successful for", eff, "effective samples.\n" ) )
  cat( paste("Dependence structure bootstrap mean parameter estimates:\n" ))
  print( co[[ 1 ]] )
  cat("\n")

	}

