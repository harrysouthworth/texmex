revTransform <-
function (x, data, qu, th = 0, sigma = 1, xi = 0, method = "mixture") {
   if (!is.element(method, c("mixture", "empirical")))
       stop("method should be 'mixture' or 'empirical'")

   n <- length(data)
   probs <- (1:n)/(n + 1)
   px <- vapply(x,
                function(x, p) {
                  p[[which.min(abs(x-p))]]
                }, 0, p=probs)
   px <- as.integer(round(px * (1 + n)))
   res <- sort(data)[px]
   
   if (method == "mixture"){
     # Real data contain ties which can cause x[res > th] < qu, res[res < th] > qu
     i.x <- x >= qu
     i.r <- res > th
     i.rx <- apply(cbind(i.x, i.r), 1, all)
     if (sum(i.rx > 0)){
       wh <- u2gpd(x[i.rx], p=1-qu, th=th, sigma=sigma, xi=xi)
       rth <- res[i.rx]
       o <- order(rth)
       rth <- rth[o]
       rth[length(rth):(length(rth) - length(wh) + 1)] <- rev(sort(wh))
       rth <- rth[order(o)]
       res[i.rx] <- rth
     }
   }

   res[order(x)] <- sort(res)
   res
}
