`mrlPlot` <-
function (data, umin = min(data), umax = max(data) - 0.1, conf = 0.95, 
    nint = 100, ...) 
{
    data <- c(data)
    x <- xu <- xl <- numeric(nint)
    u <- seq(umin, umax, length = nint)
    for (i in 1:nint) {
        data <- data[data > u[i]]
        x[i] <- mean(data - u[i])
        sdev <- sqrt(var(data))
        n <- length(data)
        xu[i] <- x[i] + (qnorm((1 + conf)/2) * sdev)/sqrt(n)
        xl[i] <- x[i] - (qnorm((1 + conf)/2) * sdev)/sqrt(n)
    }
    plot(u, x, type = "l", ylab = "Mean Excess", 
        ylim = c(min(xl[!is.na(xl)]), max(xu[!is.na(xu)])), ...)
    lines(u[!is.na(xl)], xl[!is.na(xl)], lty = 2)
    lines(u[!is.na(xu)], xu[!is.na(xu)], lty = 2)
}

