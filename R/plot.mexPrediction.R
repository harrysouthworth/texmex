`plot.predict.mex` <-
function( x, pch=c( 1, 3 ), col=c( 2, 8 ), cex=c( 1, 1 ), ask = TRUE, ... ){

    if ( is.R() ) {
      d <- dim( x$data$simulated )[[ 2 ]] -1
      if ( prod( par( "mfrow" ) ) < d ){
        if ( ask ) {
          op <- par(ask = TRUE)
          on.exit(par(op))
        }
      }
    }

	xdat <- x$data$real[, 1 ]
  upts <- seq(from =0.001,to=1-0.0001,len=100)
  xpts <- revTransform(upts,data=x$data$real[, 1 ], qu = mean(x$data$real[,1] < x$mth[1]), th=x$mth[1],sigma = x$gpd.coef[3,1], xi = x$gpd.coef[4,1])

	for( i in 2:( dim( x$data$real )[[ 2 ]] ) ){
		ydat <- x$data$real[, i ]
		xlimits <- range( xdat , x$data$simulated[ , 1 ] )
		ylimits <- range( ydat , x$data$simulated[ , i ] )

		plot( xdat , ydat , xlim=xlimits , ylim=ylimits,
			  xlab = names( x$data$simulated )[ 1 ],
			  ylab = names( x$data$simulated )[ i ],
			  type = "n",...
			 )
		points( x$data$simulated[, 1 ] , x$data$simulated[, i ] ,
				col=col[ 2 ] , pch=pch[ 2 ], cex=cex[ 2 ] )
		points( xdat, ydat , pch=pch[ 1 ], col=col[ 1 ], cex= cex[ 1 ] )
		abline( v = x$data$pth, lty=2, col=3 )
    ypts <- revTransform(upts,data=x$data$real[, i ], qu = mean(x$data$real[,i] < x$mth[i]), th=x$mth[i],sigma = x$gpd.coef[3,i], xi = x$gpd.coef[4,i])
    lines(xpts,ypts,col=3)
	}

	invisible()
}

test.plot.predict.mex <- function(){
# check reproduce Figure 6 in Heffernan and Tawn
  w <- mex(winter,mqu=0.7,penalty="none", which="NO", dqu=.7, margins="gumbel", constrain=FALSE)
  noMod <- bootmex(w)
  noPred <- predict(noMod)
  par(mfcol=c(2,2))
  res <- plot(noPred,main="Fig. 6 Heffernan and Tawn (2004)")
  checkEquals(res, NULL, msg="plot.predict.mex: correct execution")

# check for 2-d data
  R <- 20
  nsim <- 100
  wavesurge.mex <- mex(wavesurge,mqu=.7,which=1, margins="laplace")
  wavesurge.boot <- bootmex(wavesurge.mex,R=R)
  wavesurge.pred <- predict(wavesurge.boot,nsim=nsim)
  par(mfrow=c(1,1))
  res <- plot(wavesurge.pred)
  checkEquals(res, NULL, msg="plot.predict.mex: correct execution")
}
