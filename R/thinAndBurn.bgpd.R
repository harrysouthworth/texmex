thinAndBurn <- function (object, ...)
  UseMethod("thinAndBurn")

thinAndBurn.bgpd <- function(obj,burn,thin){

  if(missing(burn)){
    burn <- obj$burn
  } else {
    obj$burn <- burn
  }
  if(missing(thin)){
    thin <- obj$thin
  } else {
    obj$thin <- thin
  }
  if(is.null(obj$thin)){
    stop("thin or its reciprocal must be a positive integer, for no thinning use thin=1")
  }     
  if(is.null(obj$burn)){
    stop("burn must be a non-negative integer, for no burn in use burn=0")
  }     
  
 	if ( thin < 1 ) thin <- 1 / thin
	if ( thin %% 1 > 10^(-6) ) stop("thin, or its reciprocal, should be an integer" )
	if ( burn > dim(obj$chains)[1] ) stop( "burn-in is longer that the whole chain" )

  if (burn > 0){
     obj$param <- obj$chains[ -( 1:burn ) , ] # Remove burn-in
  } else {
     obj$param <- obj$chains
  }
	wh <- 1:dim( obj$param )[[ 1 ]] %% thin == 0
	obj$param <- obj$param[ wh , ]

  invisible(obj)
}

test(thinAndBurn.bgpd) <- function(){
require(svUnit)

# generate data to use for checking
  d <- sample(3:10,1)
  nrow <- 100
  x <- list(chains = apply(matrix(rep(1:d,each=nrow),ncol=d),2, function( o ) o*1:nrow))
  oldClass( x ) <- "bgpd" 

# test appropriate errors for misspecification of thin and burn

  checkException(thinAndBurn(x,burn=2))
  checkException(thinAndBurn(x,thin=1))
  
#  test burn in  
  burn <- sample(nrow/2,1)
  burnOnly <- thinAndBurn(x,burn=burn,thin=1)
  checkEqualsNumeric(x$chains[burn+1,], burnOnly$param[1,])

# test thinning
  thin <- 2
  thinOnly <- thinAndBurn(x,thin=thin,burn=0)

  checkEqualsNumeric(seq(thin,nrow,by=thin), thinOnly$param[,1])
  
# test thinning and burning simultaneously

  thinBurn <- thinAndBurn(x,thin=thin,burn=burn)
  
  checkEqualsNumeric(seq(burn + thin, nrow,by=thin), thinBurn$param[,1])
  
# test returned values of thin and burn

  checkEqualsNumeric(thin, thinBurn$thin,burn=0)
  checkEqualsNumeric(burn, thinBurn$burn,thin=1)
  
# test passing thin and burn via object

  x$thin <- thin
  x$burn <- burn
  thinBurn1 <- thinAndBurn(x)
  
  checkEqualsNumeric(seq(burn + thin, nrow,by=thin), thinBurn1$param[,1])
  checkEqualsNumeric(dim(thinBurn$param),dim(thinBurn1$param))
  
# test thinning and burning a previously thinned and burned object

  thin2 <- 4
  burn2 <- 4
  thinBurn2 <- thinAndBurn(thinBurn1,thin=thin2,burn=burn2)
  checkEqualsNumeric(dim(thinBurn1$chains)[2], dim(thinBurn2$chains)[2])
  checkEqualsNumeric((nrow - burn2) / thin2, dim(thinBurn2$param)[1])
}