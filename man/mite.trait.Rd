\name{mite.trait}
\alias{mite.trait}
\docType{data}
\title{Phenotype data for mites}
\description{A dataframe of sexual system and chromosome number data for mites. Contains three columns:
\describe{
  \item{Column 1}{species names matching the tip labels of \code{\link{trees.mite}}}
  \item{Column 2}{diploid chromosome numbers (continuous trait)}
  \item{Column 3}{sex determination system coded as 0 (diplodiploidy) or 1 (haplodiploidy) (discrete binary trait)}
}
Used in examples for the \code{\link{AncCond}} function.}
\examples{
data(mite.trait)
head(mite.trait)
}
\keyword{datasets}
