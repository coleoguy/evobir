\name{Normalize}
\alias{Normalize}
\title{Normalizes data between 0 and 1}
\description{Supplied with a vector, list or matrix of numbers it performs normalization.}
\usage{Normalize(x)}
\arguments{
  \item{x}{A vector list or matrix of numbers to be normalized.}
  \item{MARGIN}{indicates the nature of the normalization options are 1, 2, and 3. 1 indicates that each row should be normalized. 2 indicates that columns should be normalized. 3 indicates that the whole matrix should be normalized.}}
\details{Returns the same type of object as x with all the data normalized using the function: 
  \code{x = ((x - min(x)) / (max(x) - min(x)))}}
\references{
\url{http://coleoguy.github.io/}}
\author{Heath Blackmon}
\examples{Normalize(c(-32, 3, 100))}
\keyword{ data transformation }
\keyword{ Normalize }
