migpdCoefs <-
  # Coerce coefficients of a migpd object.
  # Useful when coefficients are coming from
  # a model with covariates and you want to
  # learn something about the dependence between
  # margins.
function(object, which, coefs){
  if (class(object) != "migpd"){
    stop("object must be of class \'migpd\'")
  }
  
  if (length(which) != length(coefs)){
    stop("which and coefs should have the same length")
  }
  
  if (length(which) == 1){
    object$models[[which]]$coefficients <- coefs
  }
  
  for (i in 1:length(which)){
    object$models[[i]]$coefficients <- coefs[[i]]
  }
  object
}

test.migpdCoefs <- function(){

  liver <- liver
  liver$ndose <- as.numeric(liver$dose)
  require(MASS,quiet=TRUE) # For rlm

  ralt <- resid(rlm(log(ALT_M) ~ log(ALT_B) + ndose, data=liver))
  rast <- resid(rlm(log(AST_M) ~ log(AST_B) + ndose, data=liver))
  ralp <- resid(rlm(log(ALP_M) ~ log(ALP_B) + ndose, data=liver))
  rtbl <- resid(rlm(log(TBL_M) ~ log(TBL_B) + ndose, data=liver))

  rliver <- data.frame(alt=ralt, ast=rast, alp=ralp, tbl=rtbl, ndose=liver$ndose)

  Dmod <- migpd(rliver[rliver$ndose == 4, 1:4], mqu=.7) # Model for dose D

  oldALTco <- coef(Dmod)[3:4, 1]

  altgpd <- gpd(alt, qu=.7, xi = ~ ndose, data=rliver)
  astgpd <- gpd(ast, qu=.7, xi = ~ ndose, data=rliver)

  altco <- c(coef(altgpd)[1], coef(altgpd)[2] + 4 * coef(altgpd)[3])
  astco <- c(coef(astgpd)[1], coef(astgpd)[2] + 4 * coef(astgpd)[3])

# Change one set of coefficients
  lmod <- migpdCoefs(Dmod, which="alt", list(altco))

  newALTco <- coef(lmod)[3:4, 1]
  newALTco[1] <- log(newALTco[1]) # For comparison with altco
  oldALTco[1] <- log(oldALTco[1])

  checkEqualsNumeric(altco, newALTco, msg="migpdCoefs: change one set of coefficients")
  checkTrue(all(newALTco != oldALTco), msg="migpdCoefs: change one set of coefficients")

# Change 2 sets of coefficients at once

  lmod <- migpdCoefs(Dmod, which=c("alt", "ast"), coef=list(altco, astco))

  newCo <- coef(lmod)[3:4, 1:2]
  oldCo <- coef(Dmod)[3:4, 1:2]

  newCo[1,] <- log(newCo[1,])
  oldCo[1,] <- log(oldCo[1,])

  checkEqualsNumeric(cbind(altco, astco),newCo, msg="migpdCoefs: change two set of coefficients at once")
  checkTrue(all(newCo != oldCo), "migpdCoefs: change two set of coefficients at once")
}
