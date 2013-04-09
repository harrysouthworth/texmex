rgpd <- function(n, sigma, xi, u = 0) {
  ## a standard GPD is (exp(xi * Exp(1)) - 1) / xi
  ## and the rest follows
  exponentials <- rexp(n)

  ## expand mu, sigma, and xi to be n long
  ## this is necessary to ensure that we get
  ## exactly n random numbers if mu, sigma, xi
  ## are greater than n long
  n     <- length(exponentials)
  sigma <- rep(sigma, length.out=n)
  xi    <- rep(xi, length.out=n)
  u     <- rep(u, length.out=n)

  u + sigma * .exprel(exponentials * xi) * exponentials
}


test.rgpd <- function(){
  ## testing an RNG...
  num.simple <- 1000
  num.quantile <- 1e6

  xi.values <- c(0, seq(-2, 2), length.out=10)
  test.quantiles <- c(0.25, 0.5, 0.75)

  core.sanity.test <- function(xi) {
    seed <- as.integer(runif(1, -1, 1)*(2**30))
    set.seed(seed)
    samples <- rgpd(num.simple, 1, xi)
    checkEquals(length(samples), num.simple,
                "rgpd: output of correct length")
    if (xi < 0) {
      checkTrue(all(samples <= -1/xi), "rgpd: upper bound check")
    }
    checkTrue(all(samples > 0), "rgpd: lower bound check")

    sigma <- rexp(1)
    mu    <- runif(1, -5, 5)
    shifted <- mu + sigma * samples
    set.seed(seed)
    checkEqualsNumeric(shifted,
                       rgpd(num.simple, sigma, xi, u=mu),
                       "rgpd: scale and shift")
  }

  quantile.test <- function(xi) {
    ## here are the sampled quantiles
    quantiles <- quantile(pgpd(rgpd(num.quantile, 1, xi),
                               1, xi),
                          probs=test.quantiles,
                          names=FALSE)
    ## this is a bit crude, but hey...
    checkEqualsNumeric(test.quantiles, quantiles,
                       tolerance=0.02,
                       "rgpd: quantile test")
  }

  lapply(xi.values, core.sanity.test)
  lapply(xi.values, quantile.test)
}

