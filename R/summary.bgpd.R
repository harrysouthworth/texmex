summary.bgpd <- function(x){
   co <- coef(x)
   se <- apply(x$param, 2, sd)
   res <- cbind(co, se)
   dimnames(res) <- list(names(co), c("Posterior mean", "SD"))
   oldClass(res) <- "summary.bgpd"
   res
}

print.summary.bgpd <- function(x){
   print(unclass(x))
   invisible()
}

