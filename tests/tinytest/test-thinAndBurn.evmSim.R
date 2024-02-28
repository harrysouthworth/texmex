# generate data to use for checking
d <- sample(3:10,1)
nrow <- 100
x <- list(chains = list(apply(matrix(rep(1:d,each=nrow),ncol=d),2, function(o) o*1:nrow)))
oldClass(x) <- "evmSim"

# test appropriate errors for misspecification of thin and burn

expect_error(thinAndBurn(x, burn=2), info="thinAndBurn.evmSim:errorsformisspecificationofthinandburn")
expect_error(thinAndBurn(x, thin=1), info="thinAndBurn.evmSim:errorsformisspecificationofthinandburn")

#  test burn in
burn <- sample(nrow/2,1)
burnOnly <- thinAndBurn(x, burn=burn,thin=1)
expect_equal(x$chains[[1]][burn+1, ], burnOnly$param[1, ], info="thinAndBurn.evmSim:burnin")

# test thinning
thin <- 2
thinOnly <- thinAndBurn(x,thin=thin,burn=0)

expect_equal(seq(thin, nrow, by=thin), thinOnly$param[, 1], info="thinAndBurn.evmSim:thinning")

# test thinning and burning simultaneously

thinBurn <- thinAndBurn(x,thin=thin,burn=burn)

expect_equal(seq(burn+thin, nrow, by=thin), thinBurn$param[, 1], info="thinAndBurn.evmSim:thinningandburningsimultaneously")

# test returned values of thin and burn

expect_equal(thin, thinBurn$thin, info="thinAndBurn.evmSim:testreturnedvalueofthin")
expect_equal(burn, thinBurn$burn, info="thinAndBurn.evmSim:testreturnedvalueofburn")

# test passing thin and burn via object

x$thin <- thin
x$burn <- burn
thinBurn1 <- thinAndBurn(x)

expect_equal(seq(burn+thin, nrow, by=thin), thinBurn1$param[, 1], info="thinAndBurn.evmSim:testpassingthinandburnviaobject")
expect_equal(dim(thinBurn$param), dim(thinBurn1$param), info="thinAndBurn.evmSim:testpassingthinandburnviaobject")

# test thinning and burning a previously thinned and burned object

thin2 <- 4
burn2 <- 4
thinBurn2 <- thinAndBurn(thinBurn1,thin=thin2,burn=burn2)
expect_equal(dim(thinBurn1$chains)[2], dim(thinBurn2$chains)[2], info="thinAndBurn.evmSim:thinningandburningapreviouslythinnedandburnedobject")
expect_equal((nrow-burn2)/thin2, dim(thinBurn2$param)[1], info="thinAndBurn.evmSim:thinningandburningapreviouslythinnedandburnedobject")
