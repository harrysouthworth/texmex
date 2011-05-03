`print.mexDependence` <-
function(x, ...){
	# print(x$call, ...)
    cat("Conditioning on ", x$conditioningVariable, " variable.\n", sep="")
	names(x$dqu) <- dimnames(x$parameters)[[2]]
	cat("\nThresholding quantiles for transformed data:\n")
	print(x$dqu, ...)
	cat("\nDependence structure parameter estimates:\n")
	if (any(abs(x$coefficients)[3:4, ] > 10^(-6))){
		print(x$coefficients, ...)
	}
	else {
		print(x$coefficients[1:2, ], ...)
	}
	invisible()
}

