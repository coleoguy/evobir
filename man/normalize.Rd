\name{normalize}
\alias{normalize}
\title{normalizes data between 0 and 1}
\description{Supplied with a vector, list or matrix of numbers it performs normalization.}
\usage{normalize(x)}
\arguments{
  \item{x}{A vector list or matrix of numbers to be normalized.}
  \item{MARGIN}{If 1 then rows are normalized and if 2 then columns are normalized}
}
\details{Returns the same type of object as x with all the data normalized using the function: \code{x = ((x - min(x)) / (max(x) - min(x))})}
\references{
\url{http://coleoguy.github.io/}}
\author{Heath Blackmon}
\examples{normalize(c(-32, 3, 100))}
\keyword{ data transformation }
\keyword{ normalize }
