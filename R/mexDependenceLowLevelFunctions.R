# check constraints on parameters under constrained Laplace estimation
ConstraintsAreSatisfied <- function(a,b,z,zpos,zneg,v){
   C1e <- a <= min(1, 1 - b*min(z)*v^(b-1), 1 - v^(b-1)*min(z) + min(zpos)/v) &
  	      a <= min(1, 1 - b*max(z)*v^(b-1), 1 - v^(b-1)*max(z) + max(zpos)/v)
				   
   C1o <- a <= 1 & 
          a > 1 - b * min(z) * v^(b-1) & 
          a > 1 - b * max(z) * v^(b-1) &
   		   (1 - 1/b)*(b*min(z))^(1/(1-b)) * (1-a)^(-b/(1 - b)) + min(zpos) > 0 &
			   (1 - 1/b)*(b*max(z))^(1/(1-b)) * (1-a)^(-b/(1 - b)) + max(zpos) > 0

	 C2e <- -a <= min(1, 1 + b*v^(b-1)*min(z), 1 + v^(b-1)*min(z) - min(zneg)/v) &
		      -a <= min(1, 1 + b*v^(b-1)*max(z), 1 + v^(b-1)*max(z) - max(zneg)/v)

   C2o <- -a <= 1 & 
          -a > 1 + b*v^(b-1)*min(z) & 
          -a > 1 + b*v^(b-1)*max(z) &
		     (1-1/b)*(-b*min(z))^(1/(1-b))*(1+a)^(-b/(1-b)) - min(zneg) > 0 &
		     (1-1/b)*(-b*max(z))^(1/(1-b))*(1+a)^(-b/(1-b)) - max(zneg) > 0

   if (any(is.na(c(C1e, C1o, C2e, C2o)))) {
       warning("Strayed into impossible area of parameter space")
       C1e <- C1o <- C2e <- C2o <- FALSE
   }
   
   (C1e | C1o) && (C2e | C2o)
}

# positive dependence Gumbel and pos or neg dependence Laplace neg likelihood function
PosGumb.Laplace.negloglik <- function(yex, ydep, a, b, m, s, constrain, v) {
    mu <- a * yex + m * yex^b
    sig <- s * yex^b
			   
    res <- sum(log(sig) + 0.5 * ((ydep - mu)/sig)^2)
    
    if (is.infinite(res)){
        if (res < 0){ 
          res <- -(10^40) 
        } else {
          res <- 10^40
        }
        warning("Infinite value of Q in mexDependence")
    } else if (constrain){
		   v <- v * max(yex)
			 zpos <- range(ydep - yex) # q0 & q1
			 z <- range((ydep - yex * a) / (yex^b)) # q0 & q1 
			 zneg <- range(ydep + yex) # q0 & q1
				   
			 if (!ConstraintsAreSatisfied(a,b,z,zpos,zneg,v)){
          res <- 10^40
       }
		}  
    
	  res
}

PosGumb.Laplace.negProfileLogLik <- function(yex, ydep, a, b, m, s, constrain, v) {

  Q <- function(yex, ydep, param, constrain, v, a, b){
    PosGumb.Laplace.negloglik(yex,ydep,a,b,m=param[1],s=param[2],constrain,v)
  }    

  o <- try(optim(par=c(m,s),fn=Q,method="L-BFGS-B", lower=c(-Inf, 10^(-10)), upper=c(Inf, Inf), a=a, b=b, yex = yex, ydep = ydep, constrain=constrain, 
                 v=v), silent=TRUE)
  
  if (class(o) == "try-error"){
			warning("Error in optim call from mexDependence")
			o <- as.list(o)
			o$par <- rep(NA, 2)
		} else if (o$convergence != 0) {
      warning("Non-convergence in mexDependence")
      o <- as.list(o)
      o$par <- rep(NA, 2)
    }
	   
  res <- list(profLik <- o$value, m=o$par[1], s=o$par[2])
  res
                 
}