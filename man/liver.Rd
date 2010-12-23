\name{liver}
\alias{liver}
\docType{data}
\title{Liver related laboratory data}
\description{
  Liver related laboratory data from a randomized, blind, parallel
  group clinical trial with 4 doses of a drug.
}
\usage{data(liver)}
\format{
  A data frame with 606 observations on the following 9 variables.
  \describe{
    \item{\code{ALP.B}}{Alkaline phosphatase at baseline. A numeric vector.}
    \item{\code{ALT.B}}{Alanine aminotransferase at baseline. A numeric vector.}
    \item{\code{AST.B}}{Aspartate aminotransferase at baseline. A numeric vector.}
    \item{\code{TBL.B}}{Total bilirubin at baseline. A numeric vector.}
    \item{\code{ALP.M}}{Alkaline phosphatase after treatment. A numeric vector.}
    \item{\code{ALT.M}}{Alanine aminotransferase after treatment. A numeric vector.}
    \item{\code{AST.M}}{Aspartate aminotransferase after treatment. A numeric vector.}
    \item{\code{TBL.M}}{Total bilirubin after treatment. A numeric vector.}
    \item{\code{dose}}{The treatment group (i.e. dose group). A factor with levels \code{A} \code{B} \code{C} \code{D}}
  }
}
\details{
   Dose A is the lowest dose, dose, B the next, C the next, and
   D the highest dose. The baseline values were taken prior to any
   treatment being received, and the clinical trial had a single
   post-baseline visit.
}
\source{
  AstraZeneca data on file.
}

\keyword{datasets}
