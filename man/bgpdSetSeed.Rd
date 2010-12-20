\name{bgpdSetSeed}
\alias{bgpdSetSeed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set the seed from a fited \code{bgpd} object.}
\description{
Set the seed from a fited \code{bgpd} object to ensure reproducibility of output.
}
\usage{
bgpdSetSeed(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{bgpd}, as returned by \code{gpd}
     using \code{method = "simulate"}.}
}
\details{
Sets the seed to the value used to fit the model.
}
\author{
Harry Southworth
}

\seealso{
\code{\link{gpd}}
}
\examples{
# mod <- gpd(rnorm(1000), qu=.7, method="simulate")
# bgpdSetSeed(mod)
}
