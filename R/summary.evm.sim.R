summary.evm.sim <- function(object, ...){
   co <- coef(object)
   se <- apply(object$param, 2, sd)
   res <- cbind(co, se)
   dimnames(res) <- list(names(co), c("Posterior mean", "SD"))
   res <- list(object$map$family, res)
   oldClass(res) <- "summary.evm.sim"
   res
}

print.summary.evm.sim <- function(x, ...){
   print(x[[1]], verbose=FALSE, ...)
   cat("\nPosterior summary:\n")
   print(unclass(x[[2]]))
   invisible()
}

show.summary.evm.sim <- print.summary.evm.sim
