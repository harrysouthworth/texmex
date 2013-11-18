inv_Laplace <-
function(p)
  {
    stopifnot( (p >= 0) & (p <= 1) )	
    -sign(p-1/2)*log(1-2*abs(p-1/2))    
  }
