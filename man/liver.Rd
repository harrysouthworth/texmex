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
    \item{\code{ALP_B}}{Alkaline phosphatase at baseline. A numeric vector.}
    \item{\code{ALT_B}}{Alanine aminotransferase at baseline. A numeric vector.}
    \item{\code{AST_B}}{Aspartate aminotransferase at baseline. A numeric vector.}
    \item{\code{TBL_B}}{Total bilirubin at baseline. A numeric vector.}
    \item{\code{ALP_M}}{Alkaline phosphatase after treatment. A numeric vector.}
    \item{\code{ALT_M}}{Alanine aminotransferase after treatment. A numeric vector.}
    \item{\code{AST_M}}{Aspartate aminotransferase after treatment. A numeric vector.}
    \item{\code{TBL_M}}{Total bilirubin after treatment. A numeric vector.}
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
