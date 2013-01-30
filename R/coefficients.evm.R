coefficients.evm <-
function(object, ...){
    res <- object$coefficients
#    if (length(res) == 2)
#    names(res) <- c("phi", "xi")
    res
}

coef.evm <- coefficients.evm
#function(object, ...) {
#  coefficients.evm(object)
#}

