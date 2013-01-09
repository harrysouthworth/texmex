# Author: Harry Southworth
# Date: 2013-1-9
# Purpose: Split out worker functions to perform boring tasks for extreme
#          value modelling functions.
#
###########################################################################

texmexMethod <-
    # Take character string passed by user and coerce to standard format
function(method){
    method <- casefold(method)
    if (method %in% c("o", "opt", "optim", "optimize", "optimise")){
        method <- "o"
    }
    else if (method %in% c("s", "sim", "simulate")){
        method <- "s"
    }
    else {
        stop("method should be either 'optimize' or 'simulate'")
    }
    method
}

texmexPrior <-
    # Take character string(s) passed by user and coerce to standard format
function(prior, penalty, method){
    prior <- casefold(prior)
    penalty <- casefold(penalty)
    if (prior != penalty){ # User provided one or both of prior or penalty
      if (prior == 'gaussian'){ # User provided penalty
          prior <- penalty
      }
      else if (prior != 'gaussian' & penalty != 'guassian'){ # User provided both
        stop('Provide neither or one of prior and penatly, not both.')
      }
      # else prior was provided and penalty wasn't, which is ok
    }

    if (method == 's' & !is.element(prior, c('gaussian', 'cauchy'))){
      stop('Only gaussian or cauchy prior can be used when simulating from posterior.')
    }
    prior
}

texmexTrace <-
    # Get tracing frequency for optimizer and for Markov chain
function(trace, method){
    if (method == "o"){
      if (!is.null(trace)){ # trace provided by user
          otrace <- trace
      }
      else {
          otrace <- 0
      }
    } # Close if (method == "o"
    else{ # method == "s"
      otrace <- 0
      if (is.null(trace)){
         trace <- 10000
      }
    }
    c(otrace, trace)
}
