test.extremalIndex <-
function(){
  tol <- 0.0001
  th <- quantile(rain,seq(0.7,0.99,len=10))
  for(i in 1:length(th)){
    texmex.ei <- extremalIndex(rain,threshold=th[i])
    Ferro.ei  <- sombrero:::.extRemes.exi.intervals(rain > th[i])
    
    Ferro.clust <- sombrero:::.extRemes.decluster.intervals(rain> th[i], Ferro.ei)
    texmex.clust <- declust(texmex.ei)
    
    Ferro.runs <-  sombrero:::.extRemes.decluster.runs(rain> th[i], 3)
    texmex.runs <- declust(rain,threshold=th[i],r=3,verbose=FALSE)
    
    checkEqualsNumeric(texmex.ei$EIintervals, Ferro.ei,
                       tolerance = tol,msg="extremalIndex: extRemes implementation")
    checkEqualsNumeric(texmex.clust$sizes, Ferro.clust$size,
                       tolerance = tol,msg="extremalIndex: declustering")
    
    checkEqualsNumeric(texmex.runs$nCluster,Ferro.runs$nc,msg="extremalIndex: runs declustering nc")
    checkEqualsNumeric(texmex.runs$sizes,Ferro.runs$size,msg="extremalIndex: runs declustering sizes")
  }
  
  # check passing data through data frames
  
  data <- data.frame(RAIN=rain[1:1000], rnorm=rnorm(1000), num=1:1000)
  extremalIndexRangeFit(RAIN, data,verbose=FALSE,nboot=20,nint=7)
  extremalIndexRangeFit(data$RAIN,verbose=FALSE,nboot=20,nint=7)
  
  data.de <- declust(RAIN,data=data,th=th[1],verb=FALSE)
  resp.de <- declust(data$RAIN,th=th[1],verb=FALSE)
  
  data.ei <- extremalIndex(RAIN,data=data,threshold=th[1])
  resp.ei <- extremalIndex(data$RAIN,threshold=th[1])
  
  checkEqualsNumeric(data.ei$EIintervals,resp.ei$EIintervals,tolerance=tol,msg="extremalIndex: using data frame to pass response")
  checkEqualsNumeric(data.de$clusters,resp.de$clusters,tolerance=tol,msg="extremalIndex: using data frame to pass numeric response to declustering")
  
  # test covariate fitting
  
  ei <- extremalIndex(SO2,data=winter,threshold=20)
  d <- declust(ei)
  evm(d,phi=~NO)
  
  checkEqualsNumeric(662.9508, AIC(evm(d,phi=~NO)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(662.8874, AIC(evm(d,phi=~NO2)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(651.8747, AIC(evm(d,phi=~O3)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(663.0015, AIC(evm(d,phi=~PM10)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(651.7874, AIC(evm(d,phi=~O3,xi=~NO)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(653.2512, AIC(evm(d,phi=~O3,xi=~NO2)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(653.6385, AIC(evm(d,phi=~O3,xi=~O3)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  checkEqualsNumeric(652.9238, AIC(evm(d,phi=~O3,xi=~PM10)),tolerance=tol, msg="extremalIndex: covariate fitting after declustering")
  
}
