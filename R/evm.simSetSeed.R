bgpdSetSeed <- function(x){
  if (oldClass(x) != "bgpd"){
      stop("This function expects an object of class \'bgpd\'.")
  }

  assign(".Random.seed", x$seed, envir=.GlobalEnv)

  invisible(x$seed)
}
