par(mfrow=c(2,1))
res <- gpdRangeFit(rain, umin=0, umax=50, nint=19)
plot(res, pch=16, main=c("Figure 4.2 of Coles (2001)",""), addNexcesses=FALSE)
expect_equal(names(res), c("th", "par","hi","lo","data"),
             label="gpdRangeFit: returned list with correct names")
