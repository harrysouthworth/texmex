prod <- texmex:::.specfun.safe.product

## simple values
x <- runif(10, -5, 5)
y <- runif(length(x), -5, 5)

expect_equal(prod(x, y), pmax(x * y, -1),
             info="safe product: simple values")

## complicated values
x <- c(0, 0, 1, 1)
y <- c(Inf, -Inf, Inf, -Inf)

res <- c(0, 0, Inf, -1)
expect_equal(prod(x, y), res,
             info="safe product: complicated values")


set.seed(123456)
x <- runif(10, -5, 5)
y <- c(0.5)

expect_equal(prod(x, y),
             pmax(x * y, -1))

x <- 0.5
y <- c(1,2,3,4)

expect_equal(prod(x, y),
             pmax(x * y, -1))
