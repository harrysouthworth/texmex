pgev <- function(q, mu, sigma, xi, lower.tail=TRUE, log.p=FALSE){
  ## first shift and scale
  q <- (q - mu) / sigma

  ## now set the lengths right
  n  <- max(length(q), length(xi))
  q  <- rep(q, length.out=n)
  xi <- rep(xi, length.out=n)

  ## this handles the limits correctly... I hope...
  xiq <- .specfun.safe.product(xi, q)

  res <- - q * .log1prel(xiq)
  pexp(exp(res), lower.tail=!lower.tail, log.p=log.p)
}


test.pgev <- function() {
  probabilities <- runif(10)

  xi.values <- c(0, seq(-5, 5, length.out=10))

  core.sanity.test <- function(xi) {
    randoms <- sort(rgev(10, 0, 1, xi))
    checkTrue(all(diff(pgev(randoms, 0, 1, xi)) >= 0),
              "pgev: ascending")
    checkTrue(all(diff(pgev(randoms, 0, 1, xi, lower.tail=FALSE)) <= 0),
              "pgev: descending")

    checkEqualsNumeric(log(pgev(randoms, 0, 1, xi)),
                       pgev(randoms, 0, 1, xi, log.p=TRUE),
                       "pgev: log.p (1)")
    checkEqualsNumeric(log(pgev(randoms, 0, 1, xi, lower.tail=FALSE)),
                       pgev(randoms, 0, 1, xi, lower.tail=FALSE,
                            log.p=TRUE),
                       "pgev: log.p (2)")

    mu <- runif(1, -5, 5)
    sigma <- rexp(1)

    checkEqualsNumeric(pgev(randoms, 0, 1, xi),
                       pgev(mu + sigma * randoms, mu, sigma, xi),
                       "pgev: shift and scale")
  }

  qgev.comparison <- function(xi) {
    quantiles <- qgev(probabilities, 0, 1, xi)
    my.probs  <- pgev(quantiles, 0, 1, xi)
    checkEqualsNumeric(my.probs, probabilities, "pgev: straight test")

    quantiles <- qgev(probabilities, 0, 1, xi, lower.tail=FALSE)
    my.probs  <- pgev(quantiles, 0, 1, xi, lower.tail=FALSE)
    checkEqualsNumeric(my.probs, probabilities, "pgev: lower tail")

    my.probs.2 <- pgev(quantiles, 0, 1, xi, lower.tail=TRUE)
    checkEqualsNumeric(probabilities + my.probs.2,
                       rep(1, length(probabilities)),
                       "pgev: tail flip")
  }

  lapply(xi.values, core.sanity.test)
  lapply(xi.values, qgev.comparison)
}
