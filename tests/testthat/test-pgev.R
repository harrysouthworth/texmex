context("pgev")

test_that("pgev behaves as it should", {
  skip_on_cran()
  skip_on_travis()
  probabilities <- runif(10)
  
  xi.values <- c(0, seq(-5, 5, length.out=10))
  
  core.sanity.test <- function(xi) {
    randoms <- sort(rgev(10, 0, 1, xi))
    expect_that(all(diff(pgev(randoms, 0, 1, xi)) >= 0), is_true(),
                label="pgev: ascending")
    expect_that(all(diff(pgev(randoms, 0,1,xi,lower.tail=FALSE)) <= 0),
                is_true(),
                label="pgev: descending")
    
    expect_that(log(pgev(randoms, 0, 1,xi)),
                equals(pgev(randoms, 0, 1, xi, log.p=TRUE)),
                label="pgev: log.p (1)")
    expect_that(log(pgev(randoms, 0, 1, xi, lower.tail=FALSE)),
                equals(pgev(randoms, 0, 1, xi, lower.tail=FALSE, log.p=TRUE)),
                label="pgev: log.p (2)")

    mu <- runif(1, -5, 5)
    sigma <- rexp(1)

    expect_that(pgev(randoms, 0, 1,xi),
                equals(pgev(mu + sigma * randoms, mu, sigma, xi)),
                label="pgev: shift and scale")
  } # Close core.sanity.tests
  
  qgev.comparison <- function(xi) {
    quantiles <- qgev(probabilities, 0, 1, xi)
    my.probs  <- pgev(quantiles, 0, 1, xi)
    expect_that(my.probs, equals(probabilities), "pgev:straighttest")
    
    quantiles <- qgev(probabilities, 0, 1, xi, lower.tail=FALSE)
    my.probs  <- pgev(quantiles, 0, 1, xi, lower.tail=FALSE)
    expect_that(my.probs, equals(probabilities), "pgev:lowertail")
    
    my.probs.2 <- pgev(quantiles, 0, 1, xi, lower.tail=TRUE)
    expect_that(probabilities+my.probs.2,
                equals(rep(1, length(probabilities))),
                label="pgev: tail flip")
  } # Close qgev.comparison
  
  lapply(xi.values, core.sanity.test)
  lapply(xi.values, qgev.comparison)
}
)
