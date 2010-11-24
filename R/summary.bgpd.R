summary.bgpd <- function(x){
   nms <- names(x$map$coefficients)
   co <- apply( x$param, 2, mean )
   se <- apply(x$param, 2, sd)
   res <- cbind(co, se)
   dimnames(res) <- list(nms, c("Posterior mean", "SD"))
   oldClass(res) <- "summary.bgpd"
   res
}

print.summary.bgpd <- function(x){
   print(unclass(x))
   invisible()
}

