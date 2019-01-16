context("qgev")

test_that("qgev behaves as it should", {

    ## get the probabilities that we'll use and sort them
    ## into ascending order for safekeeping
    set.seed(20101111)
    probabilities <- sort(runif(10))

    core.sanity.test <- function(xi) {
        base.quantiles <- qgev(probabilities, 0, 1, xi)
        ## check that the values are ascending
        expect_true(all(diff(base.quantiles)>=0), "qgev:ascendingquantiles")
        ## and check that we're descending correctly for upper tail
        bq2 <- qgev(probabilities, 0, 1, xi, lower.tail=FALSE)
        expect_true(all(diff(bq2)<=0), "qgev:descendingquantiles")
        ## does lower.tail work
        expect_equal(base.quantiles,
                     qgev(1 - probabilities, 0, 1, xi, lower.tail=FALSE),
                     label="qgev: lower.tail works correctly")
        ## does log.p work?
        expect_equal(base.quantiles,
                     qgev(log(probabilities), 0, 1, xi, log.p=TRUE),
                     label="qgev: log.p works")
        ## check shift and scale property
        sigma <- rexp(1)
        mu    <- runif(1, -5, 5)
        shifted <- mu + sigma * base.quantiles
        expect_equal(shifted, qgev(probabilities, mu, sigma, xi),
                     label="qgev: shift and scale")
    } # Close core.sanity.test

    lapply(c(0, seq(-5, 5, length.out=10)), core.sanity.test)

    ## known values
    expect_equal(-log(log(2)), qgev(0.5, 0, 1, 0),
                 label="qgev: median match at zero xi")
    xi <- seq(-2, 2, length=10)
    expect_equal(qgev(0.5, 0, 1, xi),
                 expm1(-log(log(2))*xi) / xi,
                 label="qgev: median match at nonzero xi")


    ############################################################################
    # The remaining tests compare against evd::qgev and are copied from the test
    # for qgpd.

    set.seed(20181011)

    evd.qgev <- .evd.qgev


    myTest <- function(mu, sig, xi, thresh, msg){
      myq <- sapply(1:nreps, function(i) qgev(x[,i], mu[i], sig[i], xi[i]))
      myp <- sapply(1:nreps, function(i) pgev(myq[,i], mu[i], sig[i], xi[i]))
      eq <- sapply(1:nreps, function(i) evd.qgev(x[,i], loc=mu[i], scale=sig[i], shape=xi[i]))
      expect_equal(eq, myq, label=paste(msg, "testusing.evd.qgev"))
      expect_equal(x, myp, label=paste(msg, "testusingqgev"))
    }

    #*************************************************************
    # Test exception for out of range probabilties
    expect_error(qgev(1.5, 1, 0, 2), label="qgev:exceptionforoutofrangeprob")
    expect_error(qgev(-1, 1, 0, 2), label="qgev:exceptionforoutofrangeprob")

    #*************************************************************
    # Test qgev. Note that .evd.qgev is NOT vectorized.

    nreps <- 100
    nsim <- 1000
    p <- matrix(runif(3 * nreps, -1, 1), ncol=3)
    p[, 2] <- p[, 2] + 1

    x <- matrix(runif(nreps * nsim), nrow=nsim)

    myTest(mu = p[, 1], sig=p[, 2], xi=p[, 3], msg="qgev: random xi")

    #*************************************************************
    # 6.5. Test qgev when some or all of xi == 0. Note that .evd.qgev is NOT vectorized.

    p[sample(1:nreps, nreps / 2), 3] <- 0
    myTest(mu=p[, 1], sig=p[, 2], xi = p[, 3], msg="qgev: some zero xi")
    p[, 3] <-  0
    myTest(mu=p[, 1], sig=p[, 2], xi = p[, 3],  msg="qgev: all zero xi")

    #*************************************************************
    # Test vectorization of qgev. Note that .evd.qgev is NOT vectorized.

    sig <- runif(nsim, 0, 2)
    xi <- runif(nsim)
    mu <- runif(nsim, -10, 10)

    x <- runif(nsim)

    myq <- qgev(x, mu, sig, xi)
    eq <- sapply(1:nsim, function(i) evd.qgev(x[i], loc=mu[i], scale=sig[i], shape=xi[i]))

    expect_equal(eq, myq, label="qgev:vectorisation")

    #*************************************************************
    # Test log.p argument

    lq <- qgev(log(x), mu, sig, xi, log.p=TRUE)
    expect_equal(myq, lq, label="qgev:log.p=TRUE")

    #*************************************************************
    # Test log.p argument

    LTq <- qgev(1-x, mu, sig, xi, lower.tail=FALSE)
    expect_equal(myq, LTq, label="qgev:lower.tail=FALSE")
})
