summary.bgpd <- function(object, ...){
   co <- coef(object)
   se <- apply(object$param, 2, sd)
   res <- cbind(co, se)
   dimnames(res) <- list(names(co), c("Posterior mean", "SD"))
   oldClass(res) <- "summary.bgpd"
   res
}

print.summary.bgpd <- function(x, ...){
   print(unclass(x))
   invisible()
}

show.summary.bgpd <- print.summary.bgpd