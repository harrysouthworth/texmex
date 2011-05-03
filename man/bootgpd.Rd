\name{bootgpd}
\alias{bootgpd}
\alias{print.bootgpd}
\alias{show.bootgpd}
\alias{summary.bootgpd}
\alias{plot.bootgpd}
\title{Parametric bootstrap for generalized Pareto models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootgpd(x, R = 99, trace = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class 'gpd'
}
  \item{R}{ The number of replicates to run. Defaults to \code{R = 99}.
}
  \item{trace}{
Report progress to the user every \code{trace} replicates. Defaults to \code{trace = 10}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, R=99, trace=10){
    if (class(x) != "gpd"){
        stop("x must be of class 'gpd'")
    }
    
    theCall <- match.call()

    nphi <- ncol(x$X.phi); nxi <- ncol(x$X.xi)
    phi <- coef(x)[1:nphi]
    xi <- coef(x)[(nphi+1):(nphi+nxi)]
    phi <- colSums(phi * t(x$X.phi))
    xi <- colSums(xi *t(x$X.xi))

    if ("priorParameters" \%in\% names(x$call)){
        pp <- x$priorParameters
    }
    else {
        pp <- list(rep(0, nphi+nxi), diag(rep(10^4, nphi+nxi)))
    }

    bfun <- function(i, xi, phi, X.phi, X.xi, co, pp, prior){
        if (i \%\% trace == 0){ cat("Replicate", i, "\n") }

        y <- rgpd(length(xi), xi=xi, sigma=exp(phi))
        gpd.fit(y=y,th=min(y) -1, X.phi=X.phi, X.xi=X.xi, penalty=prior, priorParameters=pp, start=co)$par
    }

    res <- t(sapply(1:R, bfun, xi=xi, phi=phi, X.phi=x$X.phi, X.xi=x$X.xi, co=coef(x), pp=pp, prior=x$penalty))

    res <- list(call=theCall, replicates=res)

    oldClass(res) <- "bootgpd"
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
