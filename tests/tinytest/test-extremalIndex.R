tol <- 0.0001
th <- quantile(rain,seq(0.7,0.99,len=10))
for(i in 1:length(th)){
  texmex.ei <- extremalIndex(rain,threshold=th[i])
  Ferro.ei  <- texmex:::.extRemes.exi.intervals(rain > th[i])

  Ferro.clust <- texmex:::.extRemes.decluster.intervals(rain> th[i], Ferro.ei)
  texmex.clust <- declust(texmex.ei)

  Ferro.runs <-  texmex:::.extRemes.decluster.runs(rain> th[i], 3)
  texmex.runs <- declust(rain, threshold=th[i], r=3, verbose=FALSE)

  expect_equal(texmex.ei$EIintervals, Ferro.ei, tolerance = tol,
               info="extremalIndex: extRemes implementation")
  expect_equal(texmex.clust$sizes, Ferro.clust$size, tolerance = tol,
               info="extremalIndex: declustering")

  expect_equal(texmex.runs$nCluster, Ferro.runs$nc, info="extremalIndex:runsdeclusteringnc")
  expect_equal(texmex.runs$sizes, Ferro.runs$size, info="extremalIndex:runsdeclusteringsizes")
}

# check passing data through data frames

data <- data.frame(RAIN=rain[1:1000], rnorm=rnorm(1000), num=1:1000)
extremalIndexRangeFit(RAIN, data,verbose=FALSE,nboot=10,nint=7)
extremalIndexRangeFit(data$RAIN,verbose=FALSE,nboot=10,nint=7)

data.de <- declust(RAIN, data=data, th=th[1], verb=FALSE)
resp.de <- declust(data$RAIN, th=th[1], verb=FALSE)

data.ei <- extremalIndex(RAIN,data=data,threshold=th[1])
resp.ei <- extremalIndex(data$RAIN,threshold=th[1])

expect_equal(data.ei$EIintervals, resp.ei$EIintervals, tolerance=tol,
             info="extremalIndex: using data frame to pass response")
expect_equal(data.de$clusters, resp.de$clusters, tolerance=tol,
             info="extremalIndex: using data frame to pass numeric response to declustering")

# test covariate fitting

ei <- extremalIndex(SO2,data=winter,threshold=20)
d <- declust(ei)
evm(d, phi = ~NO)

expect_equal(662.9508, as.numeric(AIC(evm(d, phi=~NO))[1]), tolerance=tol,
             info="extremalIndex: covariate fitting after declustering")
expect_equal(662.8874, as.numeric(AIC(evm(d, phi=~NO2))[1]), tolerance=tol,
             info="extremalIndex: covariate fitting after declustering")
expect_equal(651.8747, as.numeric(AIC(evm(d, phi=~O3))[1]), tolerance=tol,
             info="extremalIndex: covariate fitting after declustering")
expect_equal(663.0015, as.numeric(AIC(evm(d, phi=~PM10))[1]), tolerance=tol,
             info="extremalIndex: covariate fitting after declustering")
expect_equal(651.7874, as.numeric(AIC(evm(d, phi=~O3, xi=~NO))[1]), tolerance=tol,
             info="extremalIndex: covariate fitting after declustering")
expect_equal(653.2512, as.numeric(AIC(evm(d, phi=~O3, xi=~NO2))[1]), tolerance=tol,
             info="extremalIndex: covariate fitting after declustering")
expect_equal(653.6385, as.numeric(AIC(evm(d, phi=~O3, xi=~O3))[1]), tolerance=tol,
             info="extremalIndex: covariate fitting after declustering")
expect_equal(652.9238, as.numeric(AIC(evm(d, phi=~O3, xi=~PM10))[1]), tolerance=tol,
             info="extremalIndex: covariate fitting after declustering")
