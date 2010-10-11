`revGumbel` <-
function (x, data, qu, th = 0, sigma = 1, xi = 0, method = "mixture") {
   if (!is.element(method, c("mixture", "empirical")))
       stop("method should be 'mixture' or 'empirical'")
   p <- 1 - qu
   n <- length(data)
   probs <- (1:n)/(n + 1)
   px <- sapply(x, function(x, p) p[abs(x - p) == min(abs(x - p))], p = probs)
   px <- as.integer(round(px * (1 + n)))
   res <- sort(data)[px]
   if (method == "mixture") {
       hi <- u2gpd(x[res > th], p, th = th, sigma = sigma, xi = xi)
       res[res > th] <- hi
   }
   res[order(x)] <- sort(res)
   res
}

