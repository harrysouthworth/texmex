\name{winter}
\alias{winter}
\docType{data}
\title{ Air pollution data, winter months. }
\description{
 Air pollution data from Leeds (U.K.) city centre, collected from 1994
to 1998. The \code{winter} data set corresponds to the months
of November to February inclusive. Some outliers have been
removed, as discussed by Heffernan and Tawn, 2004.
}
\usage{data(winter)}
\format{
  A data frame with 532 observations on the following 5 variables.
  \describe{
    \item{\code{O3}}{Daily maximum ozone in parts per billion.}
    \item{\code{NO2}}{Daily maximum NO2 in parts per billion.}
    \item{\code{NO}}{Daily maximum NO in parts per billion.}
    \item{\code{SO2}}{Daily maximum SO2 in parts per billion.}
    \item{\code{PM10}}{Daily maximum PM10 in micrograms/metre^3}
  }
}
\source{
  Provided as online supplementary material to Heffernan and Tawn, 2004:
http://www.blackwellpublishing.com/rss/Readmefiles/heffernan.htm
}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\examples{
data(winter)
}
\keyword{datasets}
