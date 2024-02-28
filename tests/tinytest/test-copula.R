fun <- function(d) apply(d,2,function(x)(1:n)[rank(x)])/(1+n)
n <- 200

u2 <- cbind(sample(n),sample(n))
d2 <- fun(u2)

u3 <- cbind(sample(n),sample(n),sample(n))
d3 <- fun(u3)

expect_equal(d2, copula(u2)$copula, label="copula:2dimensional")
expect_equal(d3, copula(u3)$copula, label="copula:3dimensional")

expect_error(copula(TRUE))
expect_error(copula("text"))

dat <- data.frame(x=letters, stringsAsFactors=FALSE)
expect_error(copula(dat))

dat <- data.frame(x=letters,
                  y=seq_along(letters),
                  stringsAsFactors=FALSE)
expect_warning(copula(dat))
