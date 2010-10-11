\name{edf}
\alias{edf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Compute the empirical distribution function
}
\usage{
edf(x, na.last = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A numeric vector }
  \item{na.last}{ How to treat missing values. See \code{?rank}
                  for details. }
}
\value{
  A vector of quantiles relating to the observations in 
  \code{x}.
}

\author{ Harry Southworth }
\seealso{ \code{\link{copula}} }
\examples{
plot(winter$NO, edf(winter$NO))   
## The function is currently defined as
function (x, na.last = NA) 
{
    res <- rank(x)/(length(x) + 1)
    oldClass(res) <- "edf"
    invisible(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }

