`migpd` <-
function (data, th, qu, penalty = "gaussian", maxit = 10000,
   trace = 0, verbose=FALSE, priorParameters = NULL)
{
   theCall <- match.call()
   d <- dim(data)[2]
   if (missing(th) & missing(qu))
       stop("you must provide one of th or qu")
   if (!missing(th) & !missing(qu))
       stop("you must provide precisely one of th or qu")
   if (!missing(th))
       th <- rep(th, length = d)
   if (!missing(qu))
       qu <- rep(qu, length = d)
   if (missing(qu))
       qu <- sapply(1:d, function(i, x, th) 1 - mean(x[, i] >
           th[i]), x = data, th = th)
   if (missing(th))
       th <- sapply(1:d, function(i, x, prob) quantile(x[, i],
           prob = prob[i]), x = data, prob = qu)
   if (penalty %in% c("quadratic", "gaussian") & missing(priorParameters)) {
       gp = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2))
       priorParameters <- vector("list", length = length(th))
       for (i in 1:length(th)) priorParameters[[i]] <- gp
       names(priorParameters) <- dimnames(data)[[2]]
   }
   if (penalty %in% c("quadratic", "gaussian")) {
       nm <- names(priorParameters)
       if (is.null(nm))
           stop("priorParameters must be a named list")
       else if (any(!is.element(nm, dimnames(data)[[2]])))
           stop("the names of priorParameters must match the column names of the data")
   }
   
   wrapgpd <-function(i, x, th, penalty, maxit, verbose, trace, priorParameters) {
       if (verbose)
           cat("Fitting model", i, "\n")
       if (!is.null(priorParameters))
           priorParameters <- priorParameters[[(1:length(priorParameters))[names(priorParameters) ==
               dimnames(x)[[2]][i]]]]
       x <- c(x[, i])
       th <- th[i]
   
       gpd(x, th=th, penalty=penalty, priorParameters=priorParameters, maxit=maxit, trace=trace)
       }
       
   modlist <- lapply(1:d, wrapgpd, x=data, penalty=penalty, th=th, verbose=verbose,
                        priorParameters=priorParameters,maxit=maxit,trace=trace)
   if (length(dimnames(data)[[2]]) == dim(data)[[2]]){
       names(modlist) <- dimnames(data)[[2]]
       }
   names(th) <- names(qu) <- dimnames(data)[[2]]
   res <- list(call = theCall, models = modlist, data = data,
       th = th, qu = qu, penalty = penalty, priorParameters = priorParameters)
   res <- mexGumbel(res)
   oldClass(res) <- "migpd"
   invisible(res)
}

test(migpd) <- function(){

# values from Heffernan and Tawn (2004) Table 4. 
# Note values in published Table 4 for u_{Xi} in cols NO2 and NO Winter were reversed.

  htsummer <- rbind(qu=c(43, 43, 66.1, 22, 46),
    th = c(.9, .7, .7, .85, .7),
    sigma = c(15.8, 9.1, 32.2, 42.9, 22.8),
    xi = c(-.29, .01, .02, .08, .02))
  
  htwinter <- rbind(qu=c(23, 49, 151.6, 23, 53),
    th = rep(.7, 5),
    sigma = c(6.2, 9.3, 117.4, 19.7, 37.5),
    xi = c(-.37, -.03, -.09, .11, -.2))

  summer.gpd <- summary(migpd(summer, qu=htsummer[2,],penalty="none"),verbose=FALSE)
  winter.gpd <- summary(migpd(winter, qu=htwinter[2,],penalty="none"),verbose=FALSE)
  
  tol <- c(1,0.05,0.5,0.5)
  for(i in 1:4){
    checkEqualsNumeric(htsummer[i,], summer.gpd[i,],tol=tol[i])
    checkEqualsNumeric(htwinter[i,], winter.gpd[i,],tol=tol[i])
  }
}