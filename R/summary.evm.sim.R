summary.evm.sim <- function(object, ...){
   co <- coef(object)
   se <- apply(object$param, 2, sd)
   res <- cbind(co, se)
   dimnames(res) <- list(names(co), c("Posterior mean", "SD"))
   oldClass(res) <- "summary.bgpd"
   res
}

print.summary.evm.sim <- function(x, ...){
   print(unclass(x))
   invisible()
}

show.summary.evm.sim <- print.summary.evm.sim
