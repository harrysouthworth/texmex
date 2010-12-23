bgpdSetSeed <- function(x){
  if (oldClass(x) != "bgpd"){
      stop("This function expects an object of class \'bgpd\'.")
  }
  if (is.R()){
	  assign(".Random.seed", x$seed, envir=.GlobalEnv)
  }
  else {
  	  assign(".Random.seed", x$seed, where = 1)
  }
  invisible(x$seed)
}
