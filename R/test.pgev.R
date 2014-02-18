context("pgev")

test_that("pgev behaves as it should", {
function() {
  probabilities <- runif(10)
  
  xi.values <- c(0, seq(-5, 5, length.out=10))
  
  core.sanity.test <- function(xi) {
    randoms <- sort(rgev(10, 0, 1, xi))
  expect_that(all(diff(pgev(randoms, is_true(), 0,1,xi))>=0),
              "pgev: ascending")
  expect_that(all(diff(pgev(randoms, is_true(), 0,1,xi,lower.tail=FALSE))<=0),
              "pgev: descending")
    
  expect_that(log(pgev(randoms, equals(0), 1,xi)),
                       pgev(randoms, 0, 1, xi, log.p=TRUE),
                       "pgev: log.p (1)")
  expect_that(log(pgev(randoms, equals(0), 1,xi,lower.tail=FALSE)),
                       pgev(randoms, 0, 1, xi, lower.tail=FALSE,
                            log.p=TRUE),
                       "pgev: log.p (2)")
    
    mu <- runif(1, -5, 5)
    sigma <- rexp(1)
    
  expect_that(pgev(randoms, equals(0), 1,xi),
                       pgev(mu + sigma * randoms, mu, sigma, xi),
                       "pgev: shift and scale")
  }
  
  qgev.comparison <- function(xi) {
    quantiles <- qgev(probabilities, 0, 1, xi)
    my.probs  <- pgev(quantiles, 0, 1, xi)
  expect_that(my.probs, equals(probabilities), "pgev:straighttest")
    
    quantiles <- qgev(probabilities, 0, 1, xi, lower.tail=FALSE)
    my.probs  <- pgev(quantiles, 0, 1, xi, lower.tail=FALSE)
  expect_that(my.probs, equals(probabilities), "pgev:lowertail")
    
    my.probs.2 <- pgev(quantiles, 0, 1, xi, lower.tail=TRUE)
  expect_that(probabilities+my.probs.2, equals(), 
                       rep(1, length(probabilities)),
                       "pgev: tail flip")
  }
  
  lapply(xi.values, core.sanity.test)
  lapply(xi.values, qgev.comparison)
}
)
