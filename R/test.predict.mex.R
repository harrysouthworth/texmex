test.predict.mex <-
function(){
  # reproduce Table 5 in Heffernan and Tawn 2004
  smarmod <- mex(summer, mqu=c(.9, .7, .7, .85, .7), which="NO", penalty="none", dqu=.7,margins="gumbel",constrain=FALSE)
  wmarmod <- mex(winter, mqu=.7,  penalty="none", which="NO",margins="gumbel",constrain=FALSE)
  set.seed(20111010)
  NOmodWinter <- bootmex(wmarmod,trace=101)
  NOpredWinter <- predict(NOmodWinter, nsim = 500,trace=101) # matches sample size in H+T2004
  
  NOmodSummer <- bootmex(smarmod,trace=101)
  NOpredSummer <- predict(NOmodSummer, nsim = 500,trace=101)
  
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
  
  checkEqualsNumeric(Table5summer, resSummer,tolerance=tol,msg="predict.mex: Table 5 summer data")
  checkEqualsNumeric(Table5winter, resWinter,tolerance=tol,msg="predict.mex: Table 5 winter data")
  
  checkEqualsNumeric(pointEstSummer, resSummer[1,],tolerance=tol,msg="predict.mex: point est vs boot, summer data")
  checkEqualsNumeric(pointEstWinter, resWinter[1,],tolerance=tol,msg="predict.mex: point est vs boot, winter data")
  
  # check execution for 2-d data
  
  R <- 20
  nsim <- 100
  wavesurge.mex <- mex(wavesurge,mqu=.7,dqu=0.7,margins="laplace",which=1)
  wavesurge.boot <- bootmex(wavesurge.mex,R=R,trace=R+1)
  wavesurge.pred <- predict(wavesurge.boot,nsim=nsim,trace=R+1)
  
  checkEqualsNumeric(length(wavesurge.pred$replicates),R,msg="predict.mex execution for 2-d data")
  checkEqualsNumeric(dim(wavesurge.pred$replicates[[3]]),c(nsim,2))
  checkEquals(names(wavesurge.pred$replicates[[4]]),names(wavesurge),msg="predict.mex execution for 2-d data")
  
  # check predictions Laplace estimation equal to Gumbel for large samples and high threshold
  
  tol <- 0.01
  seeds <- 20:24
  set.seed(20111004)
  n <- 100000
  mqu <- c(0,0.9)
  dqu <- 0.99
  for(i in 1:5){
    x <- rgpd(n=n,sigma=1,xi=0.1)
    y <- 5 + rexp(1,5)*x + rnorm(n,0,x/max(x))
    data <- data.frame(x=x,y=y)
    
    data.gpd <- migpd(data , mqu=mqu, penalty="none")
    lap.mex <- mexDependence(data.gpd,which=1, dqu=dqu,start=c(-0.1,0.1),PlotLikDo=FALSE,v=20)
    gum.mex <- mex(data,mqu=c(0,0.9),which=1, dqu=dqu,margins="gumbel",constrain=FALSE)
    
    set.seed(seeds[i])
    lap.pred <- predict(lap.mex,nsim=10000,trace=R+1)
    set.seed(seeds[i])
    gum.pred <- predict(gum.mex,nsim=10000,trace=R+1)
    
    lap.ans <- summary(lap.pred)$ans
    gum.ans <- summary(gum.pred)$ans
    
    checkEqualsNumeric(lap.ans,gum.ans,tolerance=tol,msg=paste("predict.mex Laplace predictions equal to Gumbel, test replicate",i))
  }
}
