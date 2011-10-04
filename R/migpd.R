`migpd` <-
function (data, mth, mqu, penalty = "gaussian", maxit = 10000,
   trace = 0, verbose=FALSE, priorParameters = NULL)
{
   theCall <- match.call()
   if(is.null(colnames(data))){
    colnames(data) <- paste(rep("Column",ncol(data)),1:ncol(data),sep="")
  }
   d <- dim(data)[2]
   if (missing(mth) & missing(mqu))
       stop("you must provide one of mth or mqu")
   if (!missing(mth) & !missing(mqu))
       stop("you must provide precisely one of mth or mqu")
   if (!missing(mth))
       mth <- rep(mth, length = d)
   if (!missing(mqu))
       mqu <- rep(mqu, length = d)
   if (missing(mqu))
       mqu <- sapply(1:d, function(i, x, mth) 1 - mean(x[, i] >
           mth[i]), x = data, mth = mth)
   if (missing(mth))
       mth <- sapply(1:d, function(i, x, prob) quantile(x[, i],
           prob = prob[i]), x = data, prob = mqu)
   if (penalty %in% c("quadratic", "gaussian") & is.null(priorParameters)) {
       gp = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2))
       priorParameters <- vector("list", length = length(mth))
       for (i in 1:length(mth)) priorParameters[[i]] <- gp
       names(priorParameters) <- dimnames(data)[[2]]
   }
   else if (penalty %in% c("quadratic", "gaussian")) {
       nm <- names(priorParameters)
       if (is.null(nm))
           stop("priorParameters must be a named list")
       else if (any(!is.element(nm, dimnames(data)[[2]])))
           stop("the names of priorParameters must match the column names of the data")
   }
   
   wrapgpd <-function(i, x, mth, penalty, maxit, verbose, trace, priorParameters) {
       if (verbose)
           cat("Fitting model", i, "\n")
       if (!is.null(priorParameters))
           priorParameters <- priorParameters[[(1:length(priorParameters))[names(priorParameters) ==
               dimnames(x)[[2]][i]]]]
       x <- c(x[, i])
       mth <- mth[i]
   
       gpd(x, th=mth, penalty=penalty, priorParameters=priorParameters, maxit=maxit, trace=trace)
       }
       
   modlist <- lapply(1:d, wrapgpd, x=data, penalty=penalty, mth=mth, verbose=verbose,
                        priorParameters=priorParameters,maxit=maxit,trace=trace)
   if (length(dimnames(data)[[2]]) == dim(data)[[2]]){
       names(modlist) <- dimnames(data)[[2]]
       }
   names(mth) <- names(mqu) <- dimnames(data)[[2]]
   res <- list(call = theCall, models = modlist, data = data,
       mth = mth, mqu = mqu, penalty = penalty, priorParameters = priorParameters)
   oldClass(res) <- "migpd"
   invisible(res)
}

test.migpd <- function(){

# values from Heffernan and Tawn (2004) Table 4. 
# Note values in published Table 4 for u_{Xi} in cols NO2 and NO Winter were reversed.

  htsummer <- rbind(mqu=c(43, 43, 66.1, 22, 46),
    mth = c(.9, .7, .7, .85, .7),
    sigma = c(15.8, 9.1, 32.2, 42.9, 22.8),
    xi = c(-.29, .01, .02, .08, .02))
  
  htwinter <- rbind(mqu=c(28, 49, 151.6, 23, 53),
    mth = rep(.7, 5),
    sigma = c(6.2, 9.3, 117.4, 19.7, 37.5),
    xi = c(-.37, -.03, -.09, .11, -.2))

  summer.gpd <- summary(migpd(summer, mqu=htsummer[2,],penalty="none"),verbose=FALSE)
  winter.gpd <- summary(migpd(winter, mqu=htwinter[2,],penalty="none"),verbose=FALSE)
  
  tol <- c(1, 0.05, .5, 0.5)
  for(i in 1:4){
    checkEqualsNumeric(htsummer[i,], summer.gpd[i,],tol=tol[i],msg=paste("migpd: Table 4 summer",i))
    checkEqualsNumeric(htwinter[i,], winter.gpd[i,],tol=tol[i],msg=paste("migpd: Table 4 winter",i))
  }
  
# check excecution for 2-d data

  wavesurge.fit <- migpd(wavesurge,mq=.7)
  checkEqualsNumeric(wavesurge.fit$models$wave$loglik,gpd(wavesurge$wave,qu=0.7)$loglik,
                     tol=0.001,msg="migpd: 2-d data gpd fit wave")
  
}
