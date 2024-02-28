## so, how do we test an RNG...
num.simple <- 1000
num.quantile <- 1e6

xi.values  <- c(0, seq(-5, 5, length.out=10))
test.quantiles <- c(0.25, 0.5, 0.75)

core.sanity.test <- function(xi) {
  seed <- as.integer(runif(1, -1, 1)*(2**30))
  set.seed(seed)
  samples <- rgev(num.simple, 0, 1, xi)
  expect_equal(length(samples), num.simple,
              info = "rgev: output of correct length")
  if (xi > 0) {
    expect_true(all(samples>=-1/xi), info="rgev:lowerboundcheck")
  } else if (xi < 0) {
    expect_true(all(samples<=-1/xi), info="rgev:upperboundcheck")
  }
  ## scale and shift property
  sigma <- rexp(1)
  mu    <- runif(1, -5, 5)
  shifted <- mu + sigma * samples
  set.seed(seed)
  expect_equal(shifted, rgev(num.simple, mu, sigma, xi),info="rgev: scale and shift")
} # Close core.sanity.test

quantile.test <- function(xi) {
  ## here are the sampled quantiles
  quantiles <- quantile(pgev(rgev(num.quantile, 0, 1, xi),
                             0, 1, xi),
                        probs=test.quantiles,
                        names=FALSE)
  ## this is a bit crude, but hey...
  expect_equal(test.quantiles, quantiles, tolerance=0.02,
               scale=quantiles[2], check.attributes = F,
               info = "rgev: quantile test")
} # Close quantile.test

for (xi in xi.values){
  core.sanity.test(xi)
  quantile.test(xi)
}
