#' Joint exceedance curves
#' 
#' Calculate bivariate joint exceedance curves
#' 
#' Calculates pairs of points (x,y) for which the point exceedance probability P(X>x and Y>y) is constant.  This is available only in two dimensions: for higher dimensional data, the bivariate margin will be used and other variables ignored. Takes as input either two column matrix of observations, or output from \code{mexMonteCarlo} (in which case samples from all fitted models are used to calculate curves) or output from a call to the \code{predict} method for an object of class \code{mex} (in which case just the single fitted model is used for estimation, with the importance sample generated in the call to \code{predict} being used to calculate the joint exceedance curve).
#' @export JointExceedanceCurve
#' 
JointExceedanceCurve <- function(Sample, ExceedanceProb) {
    theCall <- match.call()
    UseMethod("JointExceedanceCurve", Sample)
}

#' @rdname JointExceedanceCurve
#' @export
JointExceedanceCurve.default <- function(Sample, ExceedanceProb,n=50,x=NULL) {
    calcJointExceedanceCurve(Sample,ExceedanceProb,n,x)
}

#' @rdname JointExceedanceCurve
#' @export
#' @method JointExceedanceCurve mexList
#' @param which Vector length two identifying which margins to use for joint exceedance curve estimation. Can be integer identifying columns of original data frame, or character vector identifying variables by name (these must match column names in original data).
JointExceedanceCurve.mexMC <- function(Sample, ExceedanceProb,n=50,x=NULL,which=1:2) {
    Sample <- as.matrix(Sample$MCsample[,which])
    Sample <- Sample[!is.na(Sample[,1]),]
    Sample <- Sample[!is.na(Sample[,2]),]

    calcJointExceedanceCurve(Sample,ExceedanceProb,n,x)
}

#' @rdname JointExceedanceCurve
#' @export
#' @method JointExceedanceCurve predict.mex
JointExceedanceCurve.predict.mex <- function(Sample, ExceedanceProb,n=50,x=NULL,which=1:2) {
    CondExceedanceProb <- 1-Sample$pqu
    Sample <- as.matrix(Sample$data$simulated[,which])
    if(ExceedanceProb > CondExceedanceProb) stop("ExceedanceProb must be less than the probability of exceeding the threshold used for importance sampling in the call to predict")
    
    calcJointExceedanceCurve(Sample,ExceedanceProb/CondExceedanceProb,n,x)
}

#' @rdname JointExceedanceCurve
#' @export calcJointExceedanceCurve
#' @param Sample Monte Carlo (or other) sample from which to calculate joint exceedance curve
#' @param ExceedanceProb Takes values between 0 and 1, constant value of joint exceedance probability for which the curve will be calculated
#' @param n If \code{x=NULL} then this is HALF the number of points at which the curve will be estimated (ie the curve is calculated at 2n locations)
#' @param x If specified by the user, the values of in the first dimension of \code{Sample} at which to calculate the curve. Defaults to \code{NULL} otherwise should be a numeric vector within the range of the first dimension of \code{Sample}.
calcJointExceedanceCurve  <- function(Sample,ExceedanceProb,n=50,x=NULL) {
    # mx, my are marginal upper limits
    # px, py are plotting points
    mx <- quantile(Sample[,1],1-ExceedanceProb)
    my <- quantile(Sample[,2],1-ExceedanceProb)
    if(is.null(x)){
        px <- seq(min(Sample[,1]),mx,length=n)
        py <- seq(min(Sample[,2]),my,length=n)
    } else {
        px <- x
    }
    
    f <- function(z,x,y) {
        sapply(z,function(z){
            g <- function(w) mean(x > z & y > w) - ExceedanceProb
            if(g(range(y)[1]) * g(range(y)[2]) < 0) {
                out <- uniroot(g,lower=range(y)[1],upper=range(y)[2])$root
            } else {
                out <- min(y[x>z])   
            }
            out
        }
        )
    }
    
    #calculate curve values at plotting points
    cx <- f(px,Sample[,1],Sample[,2])
    if(is.null(x)){
        cy <- f(py,Sample[,2],Sample[,1])
        res <- data.frame(x=c(px,cy),y=c(cx,py))
    } else {
        res <- data.frame(x=px,y=cx)
    }
    # sort results
    res <- res[order(res[,1]),]
    row.names(res) <- NULL
    oldClass(res) <- "jointExcCurve"
    res
}
