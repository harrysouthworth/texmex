\name{endPoint}
\alias{endPoint}
\alias{endPoint.gpd}
\alias{endPoint.bgpd}
\title{ Calculate upper end point for fitted GPD}
\description{
  Calculate upper end point for fitted GPD
}
\usage{
endPoint(y, verbose=TRUE,...)

\method{endPoint}{gpd}(y,verbose=TRUE,...)
\method{endPoint}{bgpd}(y,verbose=TRUE,...)

}
\arguments{
  \item{y}{ Object of class \code{gpd} or \code{bgpd}.}
  \item{verbose}{Whether to print output.}
  \item{...}{further arguments to be passed to the \code{\link{signif}} function.}
}
\value{
  In cases where the fitted shape parameter is negative, the fitted finite upper endpoint of the GPD.
}

\author{ Janet E. Heffernan}

