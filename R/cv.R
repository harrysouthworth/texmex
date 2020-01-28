#' Cross-validation for a model object
#'
#' @param object A model object.
#' @param folds The number of cross-validation folds to use. Defaults to
#'   \code{folds = 10}.
#' @param ... Other arguments to be passed through to methods.
#' @details The function is generic. At present, only objects of class 'evmOpt',
#'   as returned by \code{texmex::evm} can be used.
#' @export
cv <- function(object, folds = 10, ...){
  UseMethod("cv")
}

#' Cross-validation for the shape parameter in an extreme values model
#'
#' @param object An object of class 'evmOpt' as returned by \code{evm}.
#' @param folds Integer giving the number of cross-validation folds to use.
#'   Defaults to \code{folds = 10}.
#' @param ... Not used.
#' @param penalty String specifying the type of penalty to use. Defaults to
#'   \code{penalty = "gaussian"} which is equivalent to using a quadratic
#'   penalty. The other allowed value is \code{penalty = "lasso"} and an L1
#'   penalty is used.
#' @param range A sequence of values for the penalty parameter. Defaults to
#'   \code{range = seq(1, 25, length.out = 25)}. Must be strictly positive.
#'   The values are taken to be the reciprocals of the prior variance so mush
#'   be strictly positive.
#' @param shape String giving the name of the shape parameter. Defaults to
#'   \code{shape = NULL} and the function tries to guess.
#' @details Only the shape parameter is assumed to be penalized. The penalty
#'   can be thought of in terms of the variance of a prior distribution, which
#'   is equivalent to a quadratic penalty. Because the shape parameter will
#'   usually be between -1/2 and 1/2, a prior N(0, 1/16) distribution will
#'   likely be a good starting point, so values that span 16 will usually be
#'   appropriate.
#' @export
cv.evmOpt <- function(object, folds = 10, ..., penalty = "gaussian",
                      range = seq(1, 25, length.out = 25), shape = NULL){

  theCall <- match.call(expand.dots = TRUE)

  f <- object$family
  fn <- f$name
  parNames <- names(f$param)

  y <- object$data$y
  n <- length(y)

  if (folds > length(y)){
    stop("Number of cv folds is greater than the number of observations")
  } else if (folds < 2){
    stop("folds must be at least 2")
  }

  if (min(range) <= 0){
    stop("range must be positive")
  }

  nRange <- length(range)

  if (fn %in% c("GPD", "EGP3", "GEV")){
    shapeName <- "xi"
  } else if (fn == "CGPD"){
    shapeName <- "eta"
  } else if (is.null(shape)){
    stop("You need to tell me the name of the shape parameter")
  } else {
    shapeName <- shape
    if (!(shapeName %in% parNames))
      stop("shape not in names of family parameters")
  }

  sindex <- (1:length(parNames))[parNames == shapeName]

  dat <- data.frame(y = y, index = 1:n,
                    fold = sample(rep(1:folds, length.out = n)))

  param <- cverror <- rep(0, nRange)

  cost <- function(d, m){
    co <- matrix(coef(m), nrow = 1)
    -mean(m$family$density(d, co, model = m, log.d = TRUE))
  }

  if (penalty == "gaussian") {
    pp <- list(rep(0, length(parNames)), diag(c(rep(10^6, length(parNames) - 1), 0)))
  } else {
    pp <- list(rep(0, length(parNames)), diag(c(rep(10^(-6), length(parNames) - 1), 0)))
  }

  p <- rep(-99, times = folds)

  for (i in 1:nRange){
    pp[[2]][sindex, sindex] <- 1 / range[i]
    for (j in 1:folds){
      m <- evm(y, dat[dat$fold != j, ], family = f, th = object$threshold,
               penalty = penalty, priorParameters = pp, cov = object$cov.method)

      cverror[i] <- cverror[i] +  cost(dat$y[dat$fold == j], m)
    }
    m <- evm(y, dat, family = f, th = object$threshold,
             penalty = penalty, priorParameters = pp, cov = object$cov.method)
    param[i] <- coef(m)[paste0(shapeName, ": (Intercept)")]
  }

  res <- data.frame(penalty = range, error = cverror, estimate = param)

  res <- list(cv = res, data = dat, best = res$penalty[res$error == min(res$error)], call = theCall)
  class(res) <- "cv"
  res
}

#' @method print cv
print.cv <- function(x, ...){
  print(x$call)

  mm <- x$best
  cat(paste("\nMinimum error with penalty", round(mm, 3)))
  invisible(mm)
}

#' @method summary cv
summary.cv <- function(object, ...){
  print(object$call)

  cat("\n", max(object$data$fold), "-fold cross-validation", sep = "")
  mm <- object$best
  cat(paste("\nMinimum error with penalty", round(mm, 3)))
  invisible(mm)
}

#' @method plot cv
plot.cv <- function(x, y, ...){
  d <- x$cv
  plot(d$penalty, d$error, xlab = "Penalty", ylab = "Negative CV log-likelihood")
  invisible()
}

#' @method ggplot cv
ggplot.cv <- function(data, mapping=NULL, ..., environment = parent.frame()){
  d <- data$cv

  ggplot(d, aes(penalty, error)) +
    geom_point() +
    xlab("Penalty") + ylab("Negative CV log-likelihood")
}
