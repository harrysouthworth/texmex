# This is a weak test - it tests the structure
# of the output but not the correctness of the bootstrap coefficients; it will
# also catch ERRORs (as opposed to FAILUREs) if the code breaks.  For strong
# testing of this function, run test.predict.mex

set.seed(20120327)
smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
wmarmod <- migpd(winter, mqu=.7,  penalty="none")

mySdep <- mexDependence(smarmod, dqu=.7, which=1)
myWdep <- mexDependence(wmarmod, dqu=.7, which=1)

R <- 20

mySboot <- bootmex(mySdep, R=R,trace=R+1)
myWboot <- bootmex(myWdep, R=R,trace=R+1)

expect_equal(coef(mySdep)[[2]], mySboot$simpleDep, info="bootmex:summersimpleDepfromcallwithmarmodel")
expect_equal(coef(myWdep)[[2]], myWboot$simpleDep, info="bootmex:wintersimpleDepfromcallwithmarmodel")

expect_equal(coef(smarmod), (coef(mySboot$simpleMar)), info="bootmex:summersimpleMarfromcallwithmarmodel")
expect_equal(coef(wmarmod), (coef(myWboot$simpleMar)), info="bootmex:wintersimpleMarfromcallwithmarmodel")

expect_equal(dim(mySboot$boot[[1]]$Z)[2], (dim(mySdep$dependence$Z)[2]), info="bootmex:summerdimofresidualsfromcallwithmarmodel")
expect_equal(dim(myWboot$boot[[1]]$Z)[2], (dim(myWdep$dependence$Z)[2]), info="bootmex:winterdimofresidualsfromcallwithmarmodel")

expect_equal(dim(mySboot$boot[[1]]$dependence), (dim(coef(mySdep)[[2]])), info="bootmex:summerdimofcoefficientsfromcallwithmarmodel")
expect_equal(dim(myWboot$boot[[1]]$dependence), (dim(coef(myWdep)[[2]])), info="bootmex:winterdimofcoefficientsfromcallwithmarmodel")

expect_equal(R, (length(mySboot$boot)), info="bootmex:numberofbootstrapsamples")
expect_equal(R, (length(myWboot$boot)), info="bootmex:numberofbootstrapsamples")

expect_equal(dim(summer), (dim(mySboot$boot[[1]]$Y)), info="bootmex:sizeofbootstrapdataset")
expect_equal(dim(winter), (dim(myWboot$boot[[5]]$Y)), info="bootmex:sizeofbootstrapdataset")

smexmod <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="gumbel",constrain=FALSE, which=1)
wmexmod <- mex(winter, mqu=.7,  dqu=.7, penalty="none", margins="gumbel",constrain=FALSE, which=1)

mySboot <- bootmex(smexmod, R=R,trace=R+1)
myWboot <- bootmex(wmexmod, R=R,trace=R+1)

expect_equal(coef(smexmod)[[2]], (mySboot$simpleDep), info="bootmex:summersimpleDepfromcallwithmexmodel")
expect_equal(coef(smexmod)[[1]], (coef(mySboot$simpleMar)), info="bootmex:summersimpleMarfromcallwithmexmodel")

expect_equal(coef(wmexmod)[[2]], (myWboot$simpleDep), info="bootmex:wintersimpleDepfromcallwithmexmodel")
expect_equal(coef(wmexmod)[[1]], (coef(myWboot$simpleMar)), info="bootmex:wintersimpleMarfromcallwithmexmodel")

expect_equal(R, (length(mySboot$boot)), info="bootmex:numberofbootstrapsamples")
expect_equal(R, (length(myWboot$boot)), info="bootmex:numberofbootstrapsamples")

expect_equal(dim(summer), (dim(mySboot$boot[[1]]$Y)), info="bootmex:sizeofbootstrapdataset")
expect_equal(dim(winter), (dim(myWboot$boot[[5]]$Y)), info="bootmex:sizeofbootstrapdataset")

smexmod.1 <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="laplace",constrain=FALSE, which=1)
smexmod.2 <- mex(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none", dqu=.7, margins="laplace",constrain=TRUE,v=2, which=1)
mySboot.1 <- bootmex(smexmod.1,R=R,trace=R+1)
mySboot.2 <- bootmex(smexmod.2,R=R,trace=R+1)

expect_equal(coef(smexmod.1)[[2]], (mySboot.1$simpleDep), info="bootmex:summersimpleDepfromcallwithmexmodel,constrain=FASLE")
expect_equal(coef(smexmod.1)[[1]], (coef(mySboot.1$simpleMar)), info="bootmex:summersimpleMarfromcallwithmexmodel,constrain=FASLE")
expect_equal(coef(smexmod.2)[[2]], (mySboot.2$simpleDep), info="bootmex:summersimpleDepfromcallwithmexmodel,v=2")
expect_equal(coef(smexmod.2)[[1]], (coef(mySboot.2$simpleMar)), info="bootmex:summersimpleMarfromcallwithmexmodel,v=2")

# check execution of for 2-d data

wavesurge.fit <- migpd(wavesurge,mqu=.7)
wavesurge.mex <- mexDependence(wavesurge.fit, dqu=0.8,which=1)
R <- 20

wavesurge.boot <- bootmex(wavesurge.mex,R=R,trace=R+1)

expect_equal(dim(wavesurge.boot$boot[[1]]$Z)[2], (1), info="bootmex:executionfor2-ddata")
expect_equal(dimnames(wavesurge.boot$boot[[1]]$Z)[[2]], (names(wavesurge)[2]), info="bootmex:executionfor2-ddata")
expect_equal(length(wavesurge.boot$boot), (R), info="bootmex:executionfor2-ddata")

