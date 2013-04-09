rgev <- function(n, mu, sigma, xi){
  ## use standard GEV ~ (exp(-xi log E) - 1) / xi

  ## where E is a standard exponential

  neg.log.exp <- -log(rexp(n))

  ## expand mu, sigma, and xi to be n long
  ## this is necessary to ensure that we get
  ## exactly n random numbers if mu, sigma, xi
  ## are greater than n long
  n     <- length(neg.log.exp)
  mu    <- rep(mu, length.out=n)
  sigma <- rep(sigma, length.out=n)
  xi    <- rep(xi, length.out=n)

  ## and here we go
  standard.gev <- .exprel(neg.log.exp * xi) * neg.log.exp

  mu + sigma * standard.gev
}


test.rgev <- function() {
  ## so, how do we test an RNG...
  num.simple <- 1000
  num.quantile <- 1e6

  xi.values  <- c(0, seq(-5, 5, length.out=10))
  test.quantiles <- c(0.25, 0.5, 0.75)

  core.sanity.test <- function(xi) {
    seed <- as.integer(runif(1, -1, 1)*(2**30))
    set.seed(seed)
    samples <- rgev(num.simple, 0, 1, xi)
    checkEquals(length(samples), num.simple,
                "rgev: output of correct length")
    if (xi > 0) {
      checkTrue(all(samples >= -1/xi), "rgev: lower bound check")
    } else if (xi < 0) {
      checkTrue(all(samples <= -1/xi), "rgev: upper bound check")
    }
    ## scale and shift property
    sigma <- rexp(1)
    mu    <- runif(1, -5, 5)
    shifted <- mu + sigma * samples
    set.seed(seed)
    checkEqualsNumeric(shifted,
                       rgev(num.simple, mu, sigma, xi),
                       "rgev: scale and shift")
  }

  quantile.test <- function(xi) {
    ## here are the sampled quantiles
    quantiles <- quantile(pgev(rgev(num.quantile, 0, 1, xi),
                               0, 1, xi),
                          probs=test.quantiles,
                          names=FALSE)
    ## this is a bit crude, but hey...
    checkEqualsNumeric(test.quantiles, quantiles,
                       tolerance=0.02,
                       "rgev: quantile test")
  }
  lapply(xi.values, core.sanity.test)
  lapply(xi.values, quantile.test)
}



