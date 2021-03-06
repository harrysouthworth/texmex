% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/texmex-package.R
\docType{data}
\name{liver}
\alias{liver}
\title{Liver related laboratory data}
\format{
A data frame with 606 observations on the following 9 variables.
\describe{

\item{ALP.B}{Alkaline phosphatase at baseline. A numeric vector.}
\item{ALT.B}{Alanine aminotransferase at baseline. A numeric vector.}
\item{AST.B}{Aspartate aminotransferase at baseline. A numeric vector.}
\item{TBL.B}{Total bilirubin at baseline. A numeric vector.}
\item{ALP.M}{Alkaline phosphatase after treatment. A numeric vector.}
\item{ALT.M}{Alanine aminotransferase after treatment. A numeric vector.}
\item{AST.M}{Aspartate aminotransferase after treatment. A numeric vector.}
\item{TBL.M}{Total bilirubin after treatment. A numeric vector.}
\item{dose}{The treatment group (i.e. dose group). A factor with levels \code{A} \code{B} \code{C} \code{D}}
}
}
\source{
AstraZeneca data on file.
}
\usage{
data(liver)
}
\description{
Liver related laboratory data from a randomized, blind, parallel group
clinical trial with 4 doses of a drug.
}
\details{
Dose A is the lowest dose, dose, B the next, C the next, and D the highest
dose. The baseline values were taken prior to any treatment being received,
and the clinical trial had a single post-baseline visit.
}
\keyword{datasets}
