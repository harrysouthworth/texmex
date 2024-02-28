## testing an RNG...
num.simple <- 1000
num.quantile <- 1e6

xi.values <- c(0, seq(-2, 2), length.out=10)
test.quantiles <- c(0.25, 0.5, 0.75)

core.sanity.test <- function(xi) {
  seed <- as.integer(runif(1, -1, 1)*(2**30))
  set.seed(seed)
  samples <- rgpd(num.simple, 1, xi)
  expect_equal(length(samples), num.simple,
              info = "rgpd: output of correct length")
  if (xi < 0) {
    expect_true(all(samples<=-1/xi), info = "rgpd:upperboundcheck")
  }
  expect_true(all(samples>0), info = "rgpd:lowerboundcheck")

  sigma <- rexp(1)
  mu    <- runif(1, -5, 5)
  shifted <- mu + sigma * samples
  set.seed(seed)
  expect_equal(shifted, rgpd(num.simple, sigma, xi, u=mu),
               info = "rgpd: scale and shift")
}

quantile.test <- function(xi) {
  ## here are the sampled quantiles
  quantiles <- quantile(pgpd(rgpd(num.quantile, 1, xi),
                             1, xi),
                        probs=test.quantiles,
                        names=FALSE)
  ## this is a bit crude, but hey...
  expect_equal(test.quantiles, quantiles, tolerance=0.02,
              info = "rgpd: quantile test")
}

for (xi in xi.values){
  core.sanity.test(xi)
  quantile.test(xi)
}
