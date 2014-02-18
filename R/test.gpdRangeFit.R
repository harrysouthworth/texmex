context("gpdRangeFit")

test_that("gpdRangeFit behaves as it should", {
    par(mfrow=c(2,1))
  res <- gpdRangeFit(rain, umin=0, umax=50, nint=20)
  plot(res, pch=16, main=c("Figure 4.2 of Coles (2001)",""), addNexcesses=FALSE)
  expect_that(names(res), equals(c("th"), "par","hi","lo","data"),label=pst("gpdRangeFit:returnedlistwithcorrectnames"))
}
)
