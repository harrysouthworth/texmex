.log1mexp <- texmex:::.log1mexp

x <- seq(-2, -0.01, length.out=100)
expect_equal(.log1mexp(x), log(1-exp(x)))
expect_equal(.log1mexp(-Inf), 0)

expect_true(is.na(.log1mexp(NA)))
expect_true(is.nan(.log1mexp(0/0)))
expect_true(is.nan(.log1mexp(Inf)))
