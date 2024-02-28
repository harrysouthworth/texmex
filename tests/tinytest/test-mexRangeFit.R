which <- 2
quantiles <- seq(0.5,0.9,length=5)
wmarmod <- migpd(winter, mqu=.7,  penalty="none")
wmexmod.gum <- mexDependence(wmarmod, dqu=quantiles[1], margins="gumbel", constrain=FALSE,which=which)
wmexmod.lap <- mexDependence(wmarmod, dqu=quantiles[1], margins="laplace",v=5,which=which)

R <- 3
mrf1 <- mexRangeFit(wmarmod,quantiles = quantiles,which=which,R=R,trace=R+1,v=5)
mrf2 <- mexRangeFit(wmexmod.gum,quantiles = quantiles,R=R,trace=R+1)
mrf3 <- mexRangeFit(wmexmod.lap,quantiles = quantiles,R=R,trace=R+1)

expect_error(mexRangeFit(TRUE), info="mexRangeFit:exceptionhandle")
expect_error(mexRangeFit(5, which=1), info="mexRangeFit:exceptionhandle")

expect_equal(mrf1$ests[[1]][1:2]$margins, wmexmod.lap[1:2]$margins, info="mexmod: lap")
expect_equal(mrf1$ests[[1]][1:2]$dependence$coefficients, wmexmod.lap[1:2]$dependence$coefficients, info="mexmod: lap")

expect_equal(mrf2$ests[[1]][1:2]$margins, wmexmod.gum[1:2]$margins, info="mexmod: gum")
expect_equal(mrf2$ests[[1]][1:2]$dependence$coefficients, wmexmod.gum[1:2]$dependence$coefficients, info="mexmod: gum")

expect_equal(mrf3$ests[[1]][1:2]$margins, wmexmod.lap[1:2]$margins, info="mexmod:lap2")
expect_equal(mrf3$ests[[1]][1:2]$dependence$coefficients, wmexmod.lap[1:2]$dependence$coefficients, info="mexmod:lap2")

# now 2-d data

mqu <- .7
wavesurge.fit <- migpd(wavesurge,mqu=mqu)
m <- mexDependence(wavesurge.fit,which=1,dqu=mqu)
mrf4 <- mexRangeFit(wavesurge.fit,which=1,margins="laplace",R=R,trace=R+1)
mrf5 <- mexRangeFit(m,R=R,trace=R+1)
expect_equal(mrf4$ests[[2]][1:2], mrf5$ests[[2]][1:2], info="mrf5ests")

# test specification of starting values
R <- 5
qu <- c(0.5,0.7,0.9)
mrf6 <- mexRangeFit(wavesurge.fit,which=1,margins="laplace",constrain=TRUE, start=c(0.01,0.01),R=R,trace=R+1,quantiles = qu)
mrf7 <- mexRangeFit(wavesurge.fit,which=2,margins="laplace",constrain=TRUE, start=m,R=R,trace=R+1,quantiles = qu)

# test plotting

par(mfrow=c(2,2))
plot(mrf6,main="start=(0.01,0.01)",addNexcesses=FALSE)
plot(mrf7,main=paste("start=",signif(coef(m)$dependence[1:2],2)),addNexcesses=FALSE)


