make.mvn.prior <- texmex:::.make.mvn.prior
make.quad.prior <- texmex:::.make.quadratic.penalty
make.spd.matrix <- texmex:::.random.spd.matrix
for (count in 1:100) {
  dimn <- sample(10, size=1)
  cov.matrix <- make.spd.matrix(dimn)
  centre <- rexp(dimn)
  mvnprior <- make.mvn.prior(list(centre, cov.matrix))
  point <- rexp(dimn)
  expect_equal(mvnprior(point),
               dmvnorm(point, centre, cov.matrix, log=TRUE),
               label="efficient.closures: multivariate Gaussian prior")
  quadprior <- make.quad.prior(list(centre, cov.matrix))
  expect_equal(quadprior(point),
               mahalanobis(point, centre, cov.matrix),
               label="efficient.closures: Mahalanobis distance")
  # or "A-norm", as it's otherwise called
}
