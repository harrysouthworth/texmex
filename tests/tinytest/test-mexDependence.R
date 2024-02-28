smarmod <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
wmarmod <- migpd(winter, mqu=.7,  penalty="none")

mySdepO3 <- mexDependence(smarmod,which=1,dqu=0.7,margins="gumbel",constrain=FALSE)
myWdepO3 <- mexDependence(wmarmod,which=1,dqu=0.7,margins="gumbel",constrain=FALSE)

mySdepNO2 <- mexDependence(smarmod,which=2,dqu=0.7,margins="gumbel",constrain=FALSE)
myWdepNO2 <- mexDependence(wmarmod,which=2,dqu=0.7,margins="gumbel",constrain=FALSE)

mySdepNO <- mexDependence(smarmod,which=3,dqu=0.7,margins="gumbel",constrain=FALSE)
myWdepNO <- mexDependence(wmarmod,which=3,dqu=0.7,margins="gumbel",constrain=FALSE)

mySdepSO2 <- mexDependence(smarmod,which=4,dqu=0.7,margins="gumbel",constrain=FALSE)
myWdepSO2 <- mexDependence(wmarmod,which=4,dqu=0.7,margins="gumbel",constrain=FALSE)

mySdepPM10 <- mexDependence(smarmod,which=5,dqu=0.7,margins="gumbel",constrain=FALSE)
myWdepPM10 <- mexDependence(wmarmod,which=5,dqu=0.7,margins="gumbel",constrain=FALSE)


jhSdepO3 <- matrix(c(
  0.56627103,  0.37272912, 0.0000000, 0.0000000,
  0.22029334,  0.36865296, 0.0000000, 0.0000000,
  0.28193999, -0.26731823, 0.0000000, 0.0000000,
  0.46293139, -0.23387868, 0.0000000, 0.0000000),byrow=FALSE,nrow=4)

jhSdepNO2 <- matrix(c(
  0.49290567,  0.22236302, 0.0000000, 0.0000000,
  0.38571246,  0.34379705, 0.0000000, 0.0000000,
  0.22026515, -0.17494068, 0.0000000, 0.0000000,
  0.45455612,  0.22411795, 0.0000000, 0.0000000),byrow=FALSE,nrow=4)

jhSdepNO <- matrix(c(
  0.43149222,  0.34033851, 0.0000000, 0.0000000,
  0.49992799,  0.21878814, 0.0000000, 0.0000000,
  0.19724402,  0.23839660, 0.0000000, 0.0000000,
  0.50384850,  0.18227312, 0.0000000, 0.0000000),byrow=FALSE,nrow=4)

jhSdepSO2 <- matrix(c(
  0.24400046, -0.02162792, 0.0000000, 0.0000000,
  0.08769596, -0.14758165, 0.0000000, 0.0000000,
  0.00000000, -0.04461209, 0.6865857, 0.4201682,
  0.35364948,  0.02338747, 0.0000000, 0.0000000), byrow=FALSE,nrow=4)

jhSdepPM10 <- matrix(c(
  0.08302144,  0.16604598, 0.0000000, 0.0000000,
  0.00000000,  0.57387887, 0.0000000, 0.0000000,
  0.15208086,  0.32264497, 0.0000000, 0.0000000,
  0.00000000,  0.43255493, 0.0000000, 0.0000000), byrow=FALSE,nrow=4)

jhWdepO3 <- matrix(c(
  0.00000000,  0.008072046,  0.00000000, 0.0000000,
  0.00000000,  0.034283871,  0.00000000, 0.0000000,
  0.00000000, -0.188517544,  5.14775893, 1.0000000,
  0.00000000, -0.026874734,  0.05011460, 0.1075632),byrow=FALSE,nrow=4)

jhWdepNO2 <- matrix(c(
  0.00000000, -0.553608371, -0.06047238, 0.4967213,
  0.81920276,  0.529272235,  0.00000000, 0.0000000,
  0.32246150,  0.335335739,  0.00000000, 0.0000000,
  0.85746906,  0.085265792,  0.00000000, 0.0000000),byrow=FALSE,nrow=4)

jhWdepNO <- matrix(c(
  0.00000000, -0.504344703, -1.41890419, 0.0000000,
  0.75819233,  0.378119827,  0.00000000, 0.0000000,
  0.32199902, -0.350339706,  0.00000000, 0.0000000,
  0.73227271, -0.105822435,  0.00000000, 0.0000000),byrow=FALSE,nrow=4)

jhWdepSO2 <- matrix(c(
  0.00000000, -0.485253436, -1.27253412, 0.0000000,
  0.00000000, -0.018577905,  0.63501876, 0.3862878,
  0.00000000,  0.000000000,  0.76856266, 0.4916768,
  0.03626605, -0.316472032,  0.00000000, 0.0000000),byrow=FALSE,nrow=4)

jhWdepPM10 <- matrix(c(
  0.00000000,  0.064075145,  0.00000000, 0.0000000,
  0.86288696,  0.584629421,  0.00000000, 0.0000000,
  0.59510081,  0.569002154,  0.00000000, 0.0000000,
  0.10412199,  0.207529741,  0.00000000, 0.0000000),byrow=FALSE,nrow=4)

tol <- 0.23

expect_equal(c(jhWdepO3), c(myWdepO3$dependence$coefficients[1:4, ]), tolerance=tol, info="mexDependence:WinterO3")
expect_equal(c(jhWdepNO2), c(myWdepNO2$dependence$coefficients[1:4, ]), tolerance=tol,info="mexDependence:WinterNO2")
expect_equal(c(jhWdepNO), c(myWdepNO$dependence$coefficients[1:4, ]), tolerance=tol,info="mexDependence:WinterNO")
expect_equal(c(jhWdepSO2), c(myWdepSO2$dependence$coefficients[1:4, ]), tolerance=tol,info="mexDependence:WinterSO2")
expect_equal(c(jhWdepPM10), c(myWdepPM10$dependence$coefficients[1:4, ]), tolerance=tol,info="mexDependence:WinterPM10")

expect_equal(c(jhSdepO3), c(mySdepO3$dependence$coefficients[1:4, ]), tolerance=tol,info="mexDependence:SummerO3")
expect_equal(c(jhSdepNO2), c(mySdepNO2$dependence$coefficients[1:4, ]), tolerance=tol,info="mexDependence:SummerNO2")
expect_equal(c(jhSdepNO), c(mySdepNO$dependence$coefficients[1:4, ]), tolerance=tol,info="mexDependence:SummerNO")
expect_equal(c(jhSdepSO2), c(mySdepSO2$dependence$coefficients[1:4, ]), tolerance=tol,info="mexDependence:SummerSO2")
expect_equal(c(jhSdepPM10), c(mySdepPM10$dependence$coefficients[1:4, ]), tolerance=tol,info="mexDependence:SummerPM10")

# test marginal transformations original to Gumbel and reverse

p.x <- runif(50)
q.x <- -log(-log(p.x))
expect_equal(mySdepO3$dependence$margins[[1]],"gumbel",info="mexDependence: gumbel marginal name")
expect_equal(mySdepO3$dependence$margins$p2q(p.x), q.x, info="mexDependence: gumbel marginal transform p2q")
expect_equal(mySdepO3$dependence$margins$q2p(q.x), p.x, info="mexDependence: gumbel marginal transform q2p")

# test functionality with 2-d data

wavesurge.fit <- migpd(wavesurge,mqu=.8)
dqu <- 0.8
which <- 1
wavesurge.mex <- mexDependence(wavesurge.fit,which=which,dqu=dqu)
op <- options(warn=-1) # since want it to throw a warning but still execute correctly, which is what is tested
wavesurge.dep <- mexDependence(wavesurge.fit,which=which)
options(op)

expect_equal(wavesurge.mex[[1]], wavesurge.dep[[1]], info="mexDependence:missingdquargument")

expect_equal(wavesurge.mex[[2]]$coefficients, wavesurge.dep[[2]]$coefficients, info="mexDependence:missingdquargument")
expect_equal(wavesurge.mex[[2]]$Z, wavesurge.dep[[2]]$Z, info="mexDependence:missingdquargument")
expect_equal(wavesurge.mex[[2]]$dth, wavesurge.dep[[2]]$dth, info="mexDependence:missingdquargument")
expect_equal(wavesurge.mex[[2]]$dqu, wavesurge.dep[[2]]$dqu, info="mexDependence:missingdquargument")
expect_equal(wavesurge.mex[[2]]$which, wavesurge.dep[[2]]$which, info="mexDependence:missingdquargument")
expect_equal(wavesurge.mex[[2]]$conditioningVariable, wavesurge.dep[[2]]$conditioningVariable, info="mexDependence:missingdquargument")
expect_equal(wavesurge.mex[[2]]$loglik, wavesurge.dep[[2]]$loglik, info="mexDependence:missingdquargument")

expect_equal(wavesurge.mex[[2]]$margins[[1]], wavesurge.dep[[2]]$margins[[1]], info="mexDependence:missingdquargument")

expect_equal(wavesurge.mex[[2]]$constrain, wavesurge.dep[[2]]$constrain, info="mexDependence:missingdquargument")
expect_equal(wavesurge.mex[[2]]$v, wavesurge.dep[[2]]$v, info="mexDependence:missingdquargument")

expect_equal(dim(wavesurge.mex$dependence$Z), c(578, 1), info="mexDependence:executionfor2-ddata")
expect_equal(wavesurge.mex$dependence$dqu, dqu, info="mexDependence:executionfor2-ddata")
expect_equal(wavesurge.mex$dependence$which, which, info="mexDependence:executionfor2-ddata")

# test specification of starting values

which <- 2
dqu <- 0.8
summer.fit <- migpd(summer, mqu=c(.9, .7, .7, .85, .7), penalty="none")
summer.mex1 <- mexDependence(summer.fit,which=which,dqu=dqu)
summer.mex2 <- mexDependence(summer.fit,which=which,dqu=dqu,start=summer.mex1)
summer.mex3 <- mexDependence(summer.fit,which=which,dqu=dqu,start=summer.mex1$dependence$coefficients[1:2,])
summer.mex4 <- mexDependence(summer.fit,which=which,dqu=dqu,start=c(0.01,0.03))

tol <- 0.005
expect_equal(summer.mex1$dependence$coefficients, summer.mex2$dependence$coefficients, tolerance=tol,
             info="mexDependence:summerstartingvalues:class(start)='mex'")
expect_equal(summer.mex1$dependence$coefficients, summer.mex3$dependence$coefficients, tolerance=tol,
             info="mexDependence:summerstartingvalues:start=matrix")
expect_equal(summer.mex1$dependence$coefficients, summer.mex4$dependence$coefficients, tolerance=tol,
             info="mexDependence:summerstartingvalues:start=vector")

# test marginal transformations original to Laplace and reverse

p.x <- runif(50)
q.x <- ifelse(p.x < 0.5,log(2*p.x),-log(2*(1-p.x)))
expect_equal(summer.mex1$dependence$margins[[1]], "laplace",info="mexDependence: laplace marginal name")
expect_equal(summer.mex1$dependence$margins$p2q(p.x), q.x, info="mexDependence: laplace marginal transform p2q")
expect_equal(summer.mex1$dependence$margins$q2p(q.x), p.x, info="mexDependence: laplace marginal transform q2p")


# test laplace constrained estimation against independent coding of implementation by Yiannis Papastathopoulos (see test.functions.R file)
# do all the possible pairs generated by the winter data set.  This covers negative and positive dependence with cases on and off the boundary.

set.seed(20111010)
pairs <- expand.grid(1:5,1:5)
pairs <- as.matrix(pairs[pairs[,1] != pairs[,2], 2:1])
margins <- summer.mex1$dependence$margins # mexTransform needs list containing "laplace", p2q and q2p transforms
marTransform <- c("mixture","empirical")
Dthresh <- c(1,2,2.5)
Mquant <- 0.7
v <- 10
k <- 5

mexTransform <- texmex:::mexTransform
initial_posneg <- texmex:::initial_posneg
estimate_HT_KPT_joint_posneg_nm <- texmex:::estimate_HT_KPT_joint_posneg_nm

for(marTrans in marTransform){
  for(dth in Dthresh){
    HTSest <- KTPest <- matrix(0,ncol=2,nrow=dim(pairs)[1])

    for(i in 1:dim(pairs)[1]){
      mar.model <- migpd(winter[,pairs[i,]], mqu=Mquant,penalty="none")
      mar.trans <- mexTransform(mar.model, margins=margins, method=marTrans)
      X <- list(mar.trans$transformed)
      Dqu <- 1 - mean(mar.trans$transformed[,1] > dth)
      init <- initial_posneg(D=1, listdata=X, u=dth, v=v)
      a <- estimate_HT_KPT_joint_posneg_nm(pars=init, x=dth,listr=X,params=FALSE,k=k,v=v)
      KTPest[i,] <- a$par
      b <- mexDependence(mar.trans,which=1,dqu=Dqu,margins=margins[[1]],constrain=TRUE,v=v,
                         marTransform=marTrans,start=init,nOptim=k)
      HTSest[i,] <- coef(b)$dependence[1:2]
    }

    expect_equal(KTPest, HTSest, tolerance=tol, info=paste("mexDependence:constrainedestimationthreshold",i))
  }
}



