gpd.info <- texmex:::gpd.info

lmod <- evm(ALT.M, data=liver, qu=.5, xi=~I(240*as.numeric(dose)),
            cov="numeric")
expect_true(all(sqrt(diag(solve(gpd.info(lmod))))>0),
            info="gpd.inf: SDs positive")

# Check equality to numerical approximation in big samples
set.seed(20110923)
tol <- 10^(-3)
for (i in 1:10){
  x <- rt(10000, 10)
  junk <- evm(x, qu=.9, penalty="none", cov="numeric")
  msg <- paste("gpd.info: t", i, "equality to numerical", sep="")
  expect_equal(junk$cov, solve(gpd.info(junk)), tolerance=tol,
               info=msg)

  # check estimation when we have a penalty
  gp1 <- list(c(0, 0), diag(c(10^4, .05)))
  gp2 <- list(c(0, 0), diag(c(.1, 10^4)))
  junk1 <- evm(x, qu=.9, priorParameters = gp1, cov="numeric")
  junk2 <- evm(x, qu=.9, priorParameters = gp2, cov="numeric")
  msg1 <- paste("gpd.info: t", i, "equality to numerical, penalty on xi", sep="")
  msg2 <- paste("gpd.info: t", i, "equality to numerical, penalty on phi", sep="")
  tol <- 0.01
  expect_equal(junk1$cov, solve(gpd.info(junk1)), tolerance=tol,
               info=msg1)
  expect_equal(junk2$cov, solve(gpd.info(junk2)), tolerance=tol,
               info=msg2)

  # check estimation when we have covariates
  n <- 10000
  x <- 1/runif(n)
  data <- data.frame(x=x,y=rexp(n,exp(2 + x)))

  junk3 <- evm(y,data=data,phi =~ x,th=0)
  msg3 <- paste("gpd.info: t", i, "equality to numerical, covariates in phi", sep="")
  expect_equal(junk3$cov, solve(gpd.info(junk3)), tolerance=tol, info=msg3)

  x <- runif(n,-0.5,0.5)
  data <- data.frame(x=x,y = rgpd(n,sigma = exp(3+2*x), xi=x))

  junk4 <- evm(y,data=data,phi=~x, xi = ~ x,th=0)
  msg4 <- paste("gpd.info: t",i,"equality to numerical, covariates in phi and xi",
                sep="")
  expect_equal(junk4$cov, solve(gpd.info(junk4)), tolerance=tol, info=msg4)
}

