`mrlPlot` <-
function (data, umin = min(data), umax = max(data) - 0.1, conf = 0.95, 
    nint = 100, ...) 
{
    data <- c(data)
    AllData <- data
    x <- xu <- xl <- numeric(nint)
    Threshold <- seq(umin, umax, length = nint)
    for (i in 1:nint) {
        data <- data[data > Threshold[i]]
        x[i] <- mean(data - Threshold[i])
        sdev <- sqrt(var(data))
        n <- length(data)
        xu[i] <- x[i] + (qnorm((1 + conf)/2) * sdev)/sqrt(n)
        xl[i] <- x[i] - (qnorm((1 + conf)/2) * sdev)/sqrt(n)
    }
    plot(Threshold, x, type = "l", ylab = "Mean Excess", 
        ylim = c(min(xl[!is.na(xl)]), max(xu[!is.na(xu)])), ...)
    lines(Threshold[!is.na(xl)], xl[!is.na(xl)], lty = 2)
    lines(Threshold[!is.na(xu)], xu[!is.na(xu)], lty = 2)
    rug(AllData)
    invisible()
}

test.mrlPlot <- function(){
  par(mfrow=c(1,1))
  res <- mrlPlot(rain, main="Figure 4.1 of Coles (2001)")
  checkEquals(res,NULL,msg="mrlPlot: check execution")
}
