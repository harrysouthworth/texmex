# NOT FOR END-USERS.
# THESE ARE FUNCTIONS TAKEN FROM THE ismev AND evd PACKAGES, 
# AND FROM CODE BY YIANNIS PAPASTATHOPOULOS AND ARE USED
# SOLELY FOR TESTING texmex.
# Date evd installed: 2010-12-1
# Date ismev installed: 2010-12-1
# Date Papastathopoulos code: 2011-06-22

.evd.qgpd <-
function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{
    if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 
        1) 
        stop("`p' must contain probabilities in (0,1)")
    if (min(scale) < 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    if (lower.tail) 
        p <- 1 - p
    if (shape == 0) 
        return(loc - scale * log(p))
    else return(loc + scale * (p^(-shape) - 1)/shape)
}

.ismev.gpd.fit <-
function (xdat, threshold, npy = 365, ydat = NULL, sigl = NULL, 
    shl = NULL, siglink = identity, shlink = identity, siginit = NULL, 
    shinit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, 
    ...) 
{

	if (!is.R()){
		if (missing(shlink)){
			shlink <- function(x){ x }
		}
		if (missing(siglink)){
			siglink <- function(x){ x }
		}
	}

    z <- list()
    npsc <- length(sigl) + 1
    npsh <- length(shl) + 1
    n <- length(xdat)
    z$trans <- FALSE
    if (is.function(threshold)) 
        stop("`threshold' cannot be a function")
    u <- rep(threshold, length.out = n)
    if (length(unique(u)) > 1) 
        z$trans <- TRUE
    xdatu <- xdat[xdat > u]
    xind <- (1:n)[xdat > u]
    u <- u[xind]
    in2 <- sqrt(6 * var(xdat))/pi
    in1 <- mean(xdat, na.rm = TRUE) - 0.57722 * in2
    if (is.null(sigl)) {
        sigmat <- as.matrix(rep(1, length(xdatu)))
        if (is.null(siginit)) 
            siginit <- in2
    }
    else {
        z$trans <- TRUE
        sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])
        if (is.null(siginit)) 
            siginit <- c(in2, rep(0, length(sigl)))
    }
    if (is.null(shl)) {
        shmat <- as.matrix(rep(1, length(xdatu)))
        if (is.null(shinit)) 
            shinit <- 0.1
    }
    else {
        z$trans <- TRUE
        shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])
        if (is.null(shinit)) 
            shinit <- c(0.1, rep(0, length(shl)))
    }
    init <- c(siginit, shinit)
    z$model <- list(sigl, shl)
    z$link <- deparse(substitute(c(siglink, shlink)))
    z$threshold <- threshold
    z$nexc <- length(xdatu)
    z$data <- xdatu

	if (is.R()){
	    gpd.lik <- function(a) {
	        sc <- siglink(sigmat %*% (a[seq(1, length = npsc)]))
	        xi <- shlink(shmat %*% (a[seq(npsc + 1, length = npsh)]))
	        y <- (xdatu - u)/sc
	        y <- 1 + xi * y
	        if (min(sc) <= 0) 
	            l <- 10^6
	        else {
	            if (min(y) <= 0) 
	                l <- 10^6
	            else {
	                l <- sum(log(sc)) + sum(log(y) * (1/xi + 1))
	            }
	        }
	        l
	    }
	    x <- optim(init, gpd.lik, hessian = TRUE, method = method, 
	        control = list(maxit = maxit))
	} # Close if (is.R
	
	else {
		gpd.lik <- function(a, siglink, shlink, sigmat, shmat, npsc, npsh, xdatu, th){
	        sc <- siglink(sigmat %*% (a[seq(1, length = npsc)]))
	        xi <- shlink(shmat %*% (a[seq(npsc + 1, length = npsh)]))
	        u <- th
	        y <- (xdatu - u)/sc
	        y <- 1 + xi * y
	        if (min(sc) <= 0) 
	            l <- 10^6
	        else {
	            if (min(y) <= 0) 
	                l <- 10^6
	            else {
	                l <- sum(log(sc)) + sum(log(y) * (1/xi + 1))
	            }
	        }
	        l
		} # Close gpd.lik 
	    x <- optim(init, gpd.lik, hessian = TRUE, method = method, 
	        	   control = list(maxit = maxit, ...), siglink=siglink, shlink=shlink,
	        	   sigmat=sigmat, shmat=shmat, npsc=npsc, npsh=npsh, xdatu=xdatu, th=u)
		
	} # Close else
	

    sc <- siglink(sigmat %*% (x$par[seq(1, length = npsc)]))
    xi <- shlink(shmat %*% (x$par[seq(npsc + 1, length = npsh)]))
    z$conv <- x$convergence
    z$nllh <- x$value
    z$vals <- cbind(sc, xi, u)
    if (z$trans) {
        z$data <- -log(as.vector((1 + (xi * (xdatu - u))/sc)^(-1/xi)))
    }
    z$mle <- x$par
    z$rate <- length(xdatu)/n
    z$cov <- solve(x$hessian)
    z$se <- sqrt(diag(z$cov))
    z$n <- n
    z$npy <- npy
    z$xdata <- xdat
    if (show) {
        if (z$trans) 
            print(z[c(2, 3)])
        if (length(z[[4]]) == 1) 
            print(z[4])
        print(z[c(5, 7)])
        if (!z$conv) 
            print(z[c(8, 10, 11, 13)])
    }
    class(z) <- "gpd.fit"
    invisible(z)
}

.evd.rgpd <-
function (n, loc = 0, scale = 1, shape = 0) 
{
    if (min(scale) < 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    if (shape == 0) 
        return(loc + scale * rexp(n))
    else return(loc + scale * (runif(n)^(-shape) - 1)/shape)
}


.evd.pgpd <-
function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{
    if (min(scale) <= 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    q <- pmax(q - loc, 0)/scale
    if (shape == 0) 
        p <- 1 - exp(-q)
    else {
        p <- pmax(1 + shape * q, 0)
        p <- 1 - p^(-1/shape)
    }
    if (!lower.tail) 
        p <- 1 - p
    p
}

.evd.dgpd <-
function (x, loc = 0, scale = 1, shape = 0, log = FALSE) 
{
    if (min(scale) <= 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    d <- (x - loc)/scale
    nn <- length(d)
    scale <- rep(scale, length.out = nn)
    index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
    if (shape == 0) {
        d[index] <- log(1/scale[index]) - d[index]
        d[!index] <- -Inf
    }
    else {
        d[index] <- log(1/scale[index]) - (1/shape + 1) * log(1 + 
            shape * d[index])
        d[!index] <- -Inf
    }
    if (!log) 
        d <- exp(d)
    d
}

#*************************************************************
#*************************************************************
#*************************************************************
# Code from Yiannis Papastathopoulos
#*************************************************************
#*************************************************************
#*************************************************************

##########################################################
# **********************Unconstrained                    #
#                       Optimisation******************** #
##########################################################

Profile_likelihood_HT_unc <- function (par,listr,x,silly=-10^(40))
{
  n                <- NULL
  sig              <- NULL
  sumX             <- NULL
  temp             <- NULL
  temp2            <- NULL
  z                <- list()            
  Pl               <- silly
  X                <- vector('list',length(listr))
  Y                <- vector('list',length(listr))
  Z                <- vector('list',length(listr))
  index_alpha      <- seq(1,((2*(length(listr)) ) -1),by=2)
  index_beta       <- seq(2,((2*(length(listr)) )   ),by=2)
  alpha            <- par[index_alpha]
  beta             <- par[index_beta]    
  Z                <- vector('list',length(listr))
  Zstar            <- vector('list',length(listr))
  
   
  for(i in 1:length(listr))
    {
      temp           <- as.matrix(listr[[i]])
      X[[i]]         <- temp[,1][temp[,1]>x]
      n[i]           <- length(X[[i]])
      Y[[i]]         <- temp[,2][temp[,1]>x]
      Z[[i]]         <- (Y[[i]]  - alpha[i]*X[[i]])/(X[[i]]^beta[i])
      Zstar[[i]]     <- (Y[[i]]  - X[[i]])
      sig[i]         <- (1/n[i]) * sum ((Z[[i]]-mean(Z[[i]]))^2)
      sumX[i]        <- sum(beta[i]*log(X[[i]]))
    }
  
  if(all(alpha <= 1) & all(alpha >= -1) & all(beta < 1) ) 
    {                                       
      Pl  <- sum(((-(n/2)*log (2*pi*sig)) - sumX - (n/2)))      
    }
  if((all(alpha <= 1) ==FALSE) ||  (all(alpha >= -1)==FALSE) ||
     (all(beta < 1)==FALSE) )
    {
      Pl <- silly
    }  
  z$Pl <- Pl    
  return(z$Pl)
}
 




estimate_HT <- function(list,u,pars,params=TRUE)
                                          
  {
    res  <- optim(par=pars,Profile_likelihood_HT_unc,
                  listr=list,x=u,
                  control=list(fnscale=-1,maxit=100000))        
    ifelse(params==TRUE,return(res$par), return(res))
  }





######################################################################
# ***************end*of*unconstrained******************************* #
######################################################################




Dcond <- function(x,a,b,c,d,zi,zk)
  {
    res <- a+(x^(b-1))*zi - (x^(d-1))*zk
    return(res)
  }



###########################################################################
# ----------------------------------------------------------------------- #
# Function roots() estimates the number of stationary points of D         #
# function. (See second Theorem in draft paper for DILI)                  #
# ----------------------------------------------------------------------- #
###########################################################################

roots <- function(lev,a,c,b,d,Zj,Zk) 
  {
    #--------------------------------------------------------------------
    #Children Functions to be Used in roots().
    Dderiv <-
      function(x,alj,blj,aki,bki,Zlj,Zki)
      {
        res <- (alj-aki) + (blj*Zlj*(x^(blj-1))) - (bki*Zki*(x^(bki-1)))
        return(res)
      }

    Dderiv_01 <-
      function(x,alj,blj,aki,bki,Zlj,Zki)
      {
        res <- (alj-aki) + (blj*Zlj*((-log(x))^(blj-1))) -
          (bki*Zki*((-log(x))^(bki-1)))
        return(res)
      }
    
    is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.14)  abs(x - round(x)) < tol
    #--------------------------------------------------------------------
    
    cond_na      <- NULL
    z            <- list()
    s            <- lev
    xstar        <- lev
    xdstar       <- lev
    no_of_roots  <- 0
    
    if( (b!=0) & (d!=0) & (b!=d) & (a > c) & (Zj != 0) & (Zk != 0) )
      {
        sbase        <- (d*(d-1)*Zk)/(b*(b-1)*Zj)    
        spower       <- 1/(b-d)        
        if(sbase <= 0 & is.wholenumber((spower)/2)==FALSE) s="complex" 
        if(sbase>0 || ((is.wholenumber((spower)/2)==TRUE) & spower > 0))
          {
            s=sbase^(spower)
            cond_na <- is.na(s)

            if(cond_na==TRUE || (cond_na==FALSE & (s<=lev || s == Inf)))
              s = "complex"                        
          }

        Dprimev  <- Dderiv_01(x=exp(-lev),alj=a,blj=b,aki=c,
                              bki=d,Zlj=Zj,Zki=Zk)

        Dinf     <- Dderiv_01(x=0,alj=a,blj=b,aki=c,
                              bki=d,Zlj=Zj,Zki=Zk)
        #print(paste(s))
        if((s=="complex")==FALSE)
          {
            
            Dprimes  <- Dderiv_01(x=exp(-s),alj=a,blj=b,
                                  aki=c,bki=d,Zlj=Zj,Zki=Zk)
            if(Dprimev>0 & Dprimes>0)
              {
                no_of_roots <- 0
              }            
            if(Dprimev<0 & Dinf>0)
              {
                no_of_roots <- 1
                
                xstar     <- -log(uniroot(Dderiv_01,
                                          interval=c(exp(-lev),0),
                                          a,b,c,d,Zj,Zk)$root)          
              }
            
            if( Dprimev > 0 & Dprimes < 0 & Dinf > 0)
              {
                no_of_roots <- 2
                
                xstar     <- -log(uniroot(Dderiv_01,
                                          interval=c(exp(-lev),exp(-s)),
                                          a,b,c,d,Zj,Zk)$root)
                
                xdstar    <- -log(uniroot(Dderiv_01,
                                          interval=c(exp(-s),0),
                                          a,b,c,d,Zj,Zk)$root)           
              }
            
          }
        if((s=="complex")==TRUE)
          {
            if(Dprimev>0)  no_of_roots <- 0
            
            if(Dprimev<0 & Dinf > 0 )  
              {                
                no_of_roots <- 1              
                xstar       <- -log(uniroot(Dderiv_01,
                                            interval=c(exp(-lev),0),
                                            a,b,c,d,Zj,Zk)$root)        
              }            
          }
      }
    if(b==0  ||  d==0 || (b==d) || (a==c) ||
       (Zj == 0) || (Zk == 0) )
      {
        
        if(b==0 || Zj==0)
          {
            if((d!=0 & Zk!=0) & ((a-c) > 0))
              {
                if(((a-c)/(d*Zk))>0)
                  {
                    xstar <- ((a-c)/(d*Zk))^(1/(d-1))
                    if(xstar <= lev || xstar == Inf) xstar <- lev
                    if(xstar > lev   & xstar != Inf) no_of_roots <- 1
                  }
              }
            if((d==0 || Zk ==0)) no_of_roots <- 0
          }
        
        if(d==0 || Zk==0)
          {
            if((b!=0 & Zj!=0) & ((a-c) > 0))
              {
                if(((c-a)/(b*Zj))>0)
                  {
                    xstar <- ((c-a)/(b*Zj))^(1/(b-1))
                    if( xstar <= lev || xstar == Inf ) xstar <- lev
                    if( xstar > lev   & xstar != Inf ) no_of_roots <- 1
                  }
              }
            if((d==0 || Zk ==0)) no_of_roots <- 0            
          }
        
        if( ((a-c) == 0) & (b!=0 & d!=0) & (b!=d))
          {
            if( ((d*Zk)/(b*Zj)) > 0 )
              {
                xstar <- ((d*Zk)/(b*Zj))^(1/(b-d))
                if( xstar <= lev || xstar == Inf ) xstar <- lev
                if( xstar > lev   & xstar != Inf ) no_of_roots <- 1
              }
          }        
      }    
    z$no     <- no_of_roots
    z$s      <- s
    z$xstar  <- xstar
    z$xdstar <- xdstar       
    return(z)
  }




Profile_likelihood_cd_nm_joint_D_KT_neg <- function (par,listr,x,
                                                     Zestfun,...,v,
                                                     silly=-10^(40))
  {
    n                <- NULL
    sig              <- NULL
    sumX             <- NULL
    no_of_roots      <- NULL
    no_of_roots_star <- NULL
    temp             <- NULL
    Zq               <- NULL
    Zstarq           <- NULL
    xstar            <- NULL
    xdstar           <- NULL
    s                <- NULL
    cond_alphas      <- NULL
    cond_ord_dep     <- NULL
    cond_ord_pairs   <- NULL
    vdep             <- NULL
    z                <- list()            
    Pl               <- silly
    X                <- vector('list',length(listr))
    Y                <- vector('list',length(listr))
    Z                <- vector('list',length(listr))
    Zstar            <- vector('list',length(listr))
    index_alpha      <- seq(1,((2*(length(listr)) ) -1),by=2)
    index_beta       <- seq(2,((2*(length(listr)) )   ),by=2)
    alpha            <- par[index_alpha]
    beta             <- par[index_beta]
    xstar            <- rep(v,(length(listr)-1))
    xdstar           <- rep(v,(length(listr)-1))
    xdepstar         <- rep(vdep,length(listr))
    
    for(i in 1:length(listr))
      {
        cond_alphas[i] <- ((alpha[i]) <= 1)        
        temp       <- as.matrix(listr[[i]])
        X[[i]]     <- temp[,1][temp[,1]>x]
        vdep[i]    <- max(X[[i]])
        n[i]       <- length(X[[i]])
        Y[[i]]     <- temp[,2][temp[,1]>x]
        Z[[i]]     <- (Y[[i]] - alpha[i]*X[[i]])/(X[[i]]^beta[i])
        Zstar[[i]] <- (Y[[i]] + X[[i]])
        Zq[i]      <- Zestfun(Z[[i]],...)        
        Zstarq[i]  <- Zestfun(Zstar[[i]],...)
        sig[i]     <- (1/n[i]) * sum ((Z[[i]]-mean(Z[[i]]))^2)
        sumX[i]    <- sum(beta[i]*log(X[[i]]))
      }
    
    if(all(cond_alphas==TRUE))
      {
        for(i in 1:length(listr))
          {
            temp_roots_star      <- roots(lev=v,a=alpha[i],c=-1,
                                          b=beta[i],d=0,Zj=Zq[i],            
                                          Zk=Zstarq[i])        
            xdepstar[i]          <- temp_roots_star$xstar            
          }
      }
    
    if(all(alpha <= 1) & all(alpha >= -1) & all(beta <= 1) &
       all(cond_alphas==TRUE)) 
      {                
        for(j in 1:length(listr))
          {
            cond_ord_dep[j] <-  ( -1 <=   # (alpha[j])
                                 #mark:change v to vdep[j]
                                 (min(alpha[j],(Dcond(v,alpha[j],beta[j],-1,0,
                                               Zq[j],Zstarq[j])),#-(1e-10)),
                                      (Dcond(xdepstar[j],alpha[j],beta[j],-1,
                                             0,Zq[j],
                                             Zstarq[j])))))#-(1e-10)) )))            
          }
        
        condition <- (all(cond_ord_dep==TRUE))
        
        if( condition == TRUE )
          {
            Pl  <- sum(((-(n/2)*log (2*pi*sig)) - sumX - (n/2)))# note sig is actually sigmaSquared in the normal density
          }    
        if(condition == FALSE)
          {
            Pl  <- silly
          }          
      }
    if((all(alpha <= 1) ==FALSE) ||  (all(alpha >= -1)==FALSE) ||
       (all(beta < 1)==FALSE) || (all(cond_alphas==TRUE)==FALSE))
      {
        Pl <- silly
      }
    
     z$Pl <- Pl
#    z$Zs <- Zstarq
#    z$Zq <- Zq
      
    return(z$Pl)
  }



#########################################################################
# --------------------------------------------------------------------- #
# Function Profile_likelihood_cd_nm_joint_D_KT () returns the sum of    #
# log-likelihoods (constrained under KPT) for d=length(listr)           #
# conditional distributions (i.e. likelihood of independent             #
# conditional distributions). Each component of listr,                  #
# e.g. listr[[1]] is a matrix of dimension n[i] \times 2, with first    #
# column the conditioning variable. x is the conditioning level, v      #
# is the ordering level and Zestfun is the selected as the              #
# quantile() function with additional arguments ... If                  #
# length(listr)=1 then                                                  #
# --------------------------------------------------------------------- #
#########################################################################

Profile_likelihood_cd_nm_joint_D_KT<- function (par,listr,x,
                                                Zestfun,...,v,
                                                silly=-10^(40))
  {
    n                <- NULL
    sig              <- NULL
    sumX             <- NULL
    no_of_roots      <- NULL
    no_of_roots_star <- NULL
    temp             <- NULL
    Zq               <- NULL
    Zstarq           <- NULL
    xstar            <- NULL
    xdstar           <- NULL
    s                <- NULL
    cond_alphas      <- NULL
    cond_ord_dep     <- NULL
    cond_ord_pairs   <- NULL
    vdep             <- NULL
    z                <- list()            
    Pl               <- silly
    X                <- vector('list',length(listr))
    Y                <- vector('list',length(listr))
    Z                <- vector('list',length(listr))
    Zstar            <- vector('list',length(listr))
    index_alpha      <- seq(1,((2*(length(listr)) ) -1),by=2)
    index_beta       <- seq(2,((2*(length(listr)) )   ),by=2)
    alpha            <- par[index_alpha]
    beta             <- par[index_beta]
    xstar            <- rep(v,(length(listr)-1))
    xdstar           <- rep(v,(length(listr)-1))
    xdepstar         <- rep(vdep,length(listr))
    
    for(i in 1:length(listr))
      {
        cond_alphas[i] <- ((alpha[i]) <= 1)        
        temp           <- as.matrix(listr[[i]])
        X[[i]]         <- temp[,1][temp[,1]>x]
        vdep[i]        <- max(X[[i]])
        n[i]           <- length(X[[i]])
        Y[[i]]         <- temp[,2][temp[,1]>x]
        Z[[i]]         <- (Y[[i]] - alpha[i]*X[[i]])/(X[[i]]^beta[i])
        Zstar[[i]]     <- (Y[[i]] - X[[i]])
        Zq[i]          <- Zestfun(Z[[i]],...)        
        Zstarq[i]      <- Zestfun(Zstar[[i]],...)
        sig[i]         <- (1/n[i]) * sum ((Z[[i]]-mean(Z[[i]]))^2)
        sumX[i]        <- sum(beta[i]*log(X[[i]]))
      }
    
    if(all(cond_alphas==TRUE))
      {
        for(i in 1:length(listr))
          {
            temp_roots_star      <- roots(lev=v,a=1,c=alpha[i],
                                          b=0,d=beta[i],Zj=Zstarq[i],
                                          Zk=Zq[i])        
            xdepstar[i]          <- temp_roots_star$xstar            
          }
      }
    
    if(all(alpha <= 1) & all(alpha >= -1) & all(beta < 1) &
       all(cond_alphas==TRUE)) 
      {                
        for(j in 1:length(listr))
          {
            cond_ord_dep[j] <-  ((alpha[j]) <=
                                 #mark2:change v to vdep[j]
                                 (min(1,(Dcond(v,1,0,alpha[j],beta[j],
                                               Zstarq[j],Zq[j])),#-(1e-10)),
                                      (Dcond(xdepstar[j],1,0,alpha[j],
                                             beta[j],Zstarq[j],
                                             Zq[j])))))#-(1e-10)) )))            
          }
        
        condition <- (all(cond_ord_dep==TRUE))
        
        if( condition == TRUE )
          {
            Pl  <- sum(((-(n/2)*log (2*pi*sig)) - sumX - (n/2))) # note sig is actually sigmaSquared in the normal density
          }    
        if(condition == FALSE)
          {
            Pl  <- silly
          }          
      }
    if((all(alpha <= 1) ==FALSE) ||  (all(alpha >= -1)==FALSE) ||
       (all(beta < 1)==FALSE) || (all(cond_alphas==TRUE)==FALSE))
      {
        Pl <- silly
      }
    
     z$Pl <- Pl
#    z$Zs <- Zstarq
#    z$Zq <- Zq
      
    return(z$Pl)
  }


########################################################################
# -------------------------------------------------------------------- #
#  profile_minmax_joint_posneg_KT() combines function above to yield   #
#  all constraints for the likelihood of HT. Apart from pars all       #
#  other arguments are similar as above with different names. pars     #
#  is a vector of initial parameters.                                  #
#  ------------------------------------------------------------------- #
########################################################################


profile_minmax_joint_posneg_KT <- function(pars,listdata,u,q1=0,
                                    q2=1,...,sill=-10^(40))
  {
    loglik_min     <- NULL
    loglik_max     <- NULL
    loglik_neg_min <- NULL
    loglik_neg_max <- NULL    
    loglik         <- NULL
    
    loglik_min <- Profile_likelihood_cd_nm_joint_D_KT(par=pars,
                                                      listr=listdata,
                                                      x=u,Zestfun=quantile,
                                                      probs=q1,          
                                                      silly=sill,...)
    
    loglik_max <- Profile_likelihood_cd_nm_joint_D_KT(par=pars,
                                                      listr=listdata,
                                                      x=u,Zestfun=quantile,
                                                      probs=q2,
                                                      silly=sill,...)

    loglik_neg_min <- Profile_likelihood_cd_nm_joint_D_KT_neg(par=pars,
                                                             listr=listdata,
                                                              x=u,Zestfun=quantile,
                                                              probs=q1,
                                                              silly=sill,...)

    loglik_neg_max <- Profile_likelihood_cd_nm_joint_D_KT_neg(par=pars,
                                                              listr=listdata,
                                                              x=u,Zestfun=quantile,
                                                              probs=q2,
                                                              silly=sill,...)
    
    
    if(loglik_min == sill || loglik_max == sill ||
       loglik_neg_min == sill || loglik_neg_max == sill)
      {
        return( sill )
      }
    if(loglik_min != sill  & loglik_max != sill &
       loglik_neg_min != sill & loglik_neg_max != sill)
      {
        return( loglik_max )
      }
  }




############################################################################
# -----------------------------------------------------------------------  #
#  Function initial_posneg() is a function which is used to get            #
#  "good" initial values for argument pars of function above. It           #
#  searches over a grid of parameter values to find a point such           #
#  that the log likelihood does not fall in the constraints.               #
#  ----------------------------------------------------------------------- #
############################################################################


initial_posneg<- function(D,...)
  {    
    a    <- runif((1000*D),-1,1)
    b    <- runif((1000*D),-5,0.99)
    prop <- matrix(rbind(a,b),nrow=2*D) 
        n <- 1000
        j <- 1
        Pl <- profile_minmax_joint_posneg_KT(pars=c(a[j],b[j]),...)
    while(Pl<=-10^10)
      {
        j <- j+1
        if(j<=1000)
          {
            Pl <-  profile_minmax_joint_posneg_KT(par=c(a[j],b[j]),...)
          }
        if(j>10000)break        
      }
    if(j<=1000)
      {
        return(c(a[j],b[j]))
      }
  }


initial_posneg<- function(D,...)
  {    
    a    <- runif((1000*D),-1,1)
    b    <- runif((1000*D),-3,0.99)
    prop <- matrix(rbind(a,b),nrow=2*D) 
        n <- 1000
        j <- 1
        Pl <- profile_minmax_joint_posneg_KT(pars=prop[,j],...)
    while(Pl<=-10^10)
      {
        j <- j+1
        if(j<=1000)
          {
            Pl <-  profile_minmax_joint_posneg_KT(pars=prop[,j],...)
          }
        if(j>1000)break        
      }
    if(j<=1000)
      {
        return(prop[,j])
      }
  }

########################################################################
# -------------------------------------------------------------------- #
#  Function estimate_HT_KPT_joint_posneg_nm() estimates parameters     #
#  using KPT constraints. Arguments are same as 2nd function from      #
#  top.  k is an additional parameter that allows Nelder-Mead to be    #
#  performed more than one times.                                      #
#  ------------------------------------------------------------------- #
########################################################################

estimate_HT_KPT_joint_posneg_nm<- function(pars,x,listr,params=TRUE,...,k=3)
  {    
    temp <- NULL
    temp <- optim(par=pars,profile_minmax_joint_posneg_KT,
                  listdata=listr,u=x,...,
                  method="Nelder-Mead",
                  control=list(fnscale=-1,maxit=100000))

    tempar <- temp$par
    for(j in 1:k)
      {
        temp <- optim(par=tempar,profile_minmax_joint_posneg_KT,
                      listdata=listr,u=x,...,
                      method="Nelder-Mead",
                      control=list(fnscale=-1,maxit=100000))
        tempar <- temp$par
      }

    ifelse(params==TRUE,return(tempar), return(temp))
  }


inv_Laplace <- function(p)
  {
    stopifnot( (p >= 0) & (p <= 1) )	
    -sign(p-1/2)*log(1-2*abs(p-1/2))    
  }  
