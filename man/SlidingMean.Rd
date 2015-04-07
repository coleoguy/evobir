\name{SlidingMean}
\alias{SlidingMean}
\title{
Sliding window mean}
\description{
Calculates the mean within a sliding window of a numeric vector.  Both the step size and the window size can be set by the user.}
\usage{
SlidingMean(data, window, step)
}
\arguments{
  \item{data}{a numerical vector.}
  \item{window}{an integer setting the size of the window.}
  \item{step}{an integer setting the size of step between windows.}
}
\details{
Returns a vector of numeric values representing the mean within each window.  The length will be unequal to the original data and will be determined primarily by the step size.}
\references{
\url{http://www.uta.edu/karyodb/evobiR/}}
\author{
Heath Blackmon
}

\examples{
SlidingMean(c(1,2,1,2,10,2,1,2,1,2,3,4,5,6,2,5), 3, 1)
}
\keyword{ basic stats }
\keyword{ mode }
