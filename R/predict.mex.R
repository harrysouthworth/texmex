`predict.mex` <-
function(object, which, pqu = .99, nsim = 1000, trace=10, ...){
	theCall <- match.call()
	
    # Class can be either mex or bootmex
    theClass <- class(object)[1]
    if (! theClass %in% c("mex", "bootmex")){
        stop("object must have class 'mex' or 'bootmex'")
    }
		
	if (theClass == "bootmex" ){
        which <- object$which
        migpd <- object$simpleMar
    }
    else {
        which <- object[[2]]$which
        migpd <- object[[1]]
    }
	
	################################################################
  MakeThrowData <- function(dco,z,coxi,coxmi,data){
    y <- -log( -log( runif( nsim , min=pqu ) ) )
  	z <- as.matrix(z[ sample( 1:( dim( z )[ 1 ] ), size=nsim, replace=TRUE ) ,])
    ymi <- sapply( 1:( dim( z )[[ 2 ]] ) , makeYsubMinusI, z=z, v=dco , y=y )
  
  	ui <- exp( -exp( -y ) )
    xmi <- apply( ymi, 2, function( x ) exp( -exp( -x ) ) )

    xi <- u2gpd( ui, p = 1 - migpd$mqu[ which ], th=migpd$mth[ which ], sigma=coxi[ 1 ], xi = coxi[ 2 ] )
	
  	for( i in 1:( dim( xmi )[[ 2 ]] ) ){
		  xmi[, i ] <- revGumbel( xmi[ ,i ], as.matrix(data[,-which])[, i ],
								             th = migpd$mth[ -which ][ i ],
								             qu = migpd$mqu[ -which ][ i ],
								             sigma=coxmi[ 1,i ], xi=coxmi[ 2,i ] )
	  }
    
    sim <- data.frame( xi , xmi )
    names( sim ) <- c( names( migpd$data )[ which ],names( migpd$data )[ -which ])
    sim
  }

	################################################################
  makeYsubMinusI <- function( i, z, v , y ){
			v <- v[ , i ]
			z <- z[ , i ]
			if ( !is.na( v[ 1 ] ) ){
				if( v[ 1 ] < 10^(-5) & v[ 2 ] < 0 ){
					if( v[ 4 ] < 10^(-5 ) ) d <- 0
					else d <- v[ 4 ]
					a <- v[ 3 ] - d * log( y )
				}
				else a <- v[ 1 ] * y
			} # close if( !is.na...
			else a <- NA
			a + ( y^v[ 2 ] ) * z
		}
  
    if (theClass == "bootmex"){
    	################################################################
	    # The function lfun does most of the work
    	lfun <- function( i , bo, pqu, nsim , migpd, which ){
	    	if ( i %% trace == 0 ) cat( i, "sets done\n" )
    
            res <- MakeThrowData(dco=bo[[ i ]]$dependence,z=bo[[ i ]]$Z, coxi = bo[[i]]$GPD[,which],
                                 coxmi = as.matrix(bo[[ i ]]$GPD[,-which]),
                                 data = bo[[i]]$Y)
	    	res
    	} 
	    bootRes <- lapply( 1:length( object$boot ) , lfun ,
	        			   migpd=migpd, pqu=pqu, bo = object$boot, nsim=nsim,
	        			   which = which )
	    # bootRes contains the bootstrap simulated complete vectors X on the original 
      # scale of the data, conditional on having the _which_ component above the pqu quantile.
	} # Close if (theClass == "bootmex"
    else { bootRes <- NULL }
  
	##########################################################################
	# Get a sample using the point estimates of the parameters
	# that are suggested by the data
  
    cox <- coef(migpd)[3:4, which]
    coxmi <- as.matrix(coef(migpd)[3:4, -which])

    dall <- if (theClass == "mex") { object[[2]] }
            else { mexDependence( object$simpleMar , which=which , dqu=object$dqu ) }

    sim <- MakeThrowData(dco=dall$coefficients,z=dall$Z,coxi=cox,coxmi=coxmi,data=migpd$data)
                
    m <- 1 / ( 1 - pqu ) # Need to estimate qpu quantile
 	zeta <- 1 - migpd$mqu[ which ] # Coles, page 81
	pth <- migpd$mth[ which ] + cox[ 1 ] / cox[ 2 ] * ( ( m*zeta )^cox[ 2 ] - 1 )

	data <- list( real = data.frame( migpd$data[, which ], migpd$data[, -which] ) ,
				        simulated = sim, pth=pth)
  names(data$real)[1] <- names(migpd$data)[which] 
  
  res <- list( call = theCall , replicates = bootRes, data = data,
				       which = which, pqu = pqu,
				       mth=c( migpd$mth[ which ], migpd$mth[ -which ] ),
               gpd.coef = coef(migpd)[,c(which,(1:dim(data$real)[2])[-which])])
	
	oldClass( res ) <- "predict.mex"

	res
}

test.predict.mex <- function(){
  # reproduce Table 5 in Heffernan and Tawn 2004
  smarmod <- mex(summer, mqu=c(.9, .7, .7, .85, .7), which="NO", penalty="none", dqu=.7)
  wmarmod <- mex(winter, mqu=.7,  penalty="none", which="NO")

  NOmodWinter <- bootmex(wmarmod)
  NOpredWinter <- predict(NOmodWinter, nsim = 500) # matches sample size in H+T2004

  NOmodSummer <- bootmex(smarmod)
  NOpredSummer <- predict(NOmodSummer, nsim = 500)

  Table5winter <- rbind(c(8.3, 75.4, 569.9, 44.6, 132.3),
                      c(1.2, 4.4, 45.2, 6.7, 8.2))
  Table5summer <- rbind(c(39.6,62.2,213.5,48.5,83.7),
                        c(4.3,4.3,17.5,11.8, 7.9))
                      
  dimnames(Table5winter) <- dimnames(Table5summer) <- list(c("E(x)", "SE"),
                            c("O3", "NO2", "NO", "SO2", "PM10"))
                        
  Table5summer <- Table5summer[, c("NO", "O3", "NO2", "SO2", "PM10")]
  Table5winter <- Table5winter[, c("NO", "O3", "NO2", "SO2", "PM10")]

  resSummer <- summary(NOpredSummer)$ans[1:2,]
  resWinter <- summary(NOpredWinter)$ans[1:2,]

  pointEstSummer <- apply(NOpredSummer$data$sim,2,mean)
  pointEstWinter <- apply(NOpredWinter$data$sim,2,mean)

  tol <- 0.05
  checkEqualsNumeric(Table5summer, resSummer,tol=tol,msg="predict.mex: Table 5 summer data")
  checkEqualsNumeric(Table5winter, resWinter,tol=tol,msg="predict.mex: Table 5 winter data")
  
  checkEqualsNumeric(pointEstSummer, resSummer[1,],tol=tol,msg="predict.mex: point est vs boot, summer data")
  checkEqualsNumeric(pointEstWinter, resWinter[1,],tol=tol,msg="predict.mex: point est vs boot, winter data")

# check execution for 2-d data

  R <- 20
  nsim <- 100
  wavesurge.fit <- migpd(wavesurge,mq=.7)
  wavesurge.boot <- bootmex(wavesurge.fit,which=1,R=R)
  wavesurge.pred <- predict(wavesurge.boot,nsim=nsim)

  checkEqualsNumeric(length(wavesurge.pred$replicates),R,msg="predict.mex execution for 2-d data")
  checkEqualsNumeric(dim(wavesurge.pred$replicates[[3]]),c(nsim,2))
  checkEquals(names(wavesurge.pred$replicates[[4]]),names(wavesurge),msg="predict.mex execution for 2-d data")
  
  }
