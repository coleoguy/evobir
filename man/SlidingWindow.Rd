\name{SlidingWindow}
\alias{SlidingWindow}
\title{
Sliding window analysis}
\description{
Applies a function within a sliding window of a numeric vector or matrix. Both the step size and the window size can be set by the user. For the matrix implementation the step size and window size is constrained to be the same in both the X and Y dimensions.}
\usage{
SlidingWindow(FUN, data, window, step, strict)
}
\arguments{
  \item{FUN}{a function to be applied within each window (as a character string, e.g. \code{"mean"})}
  \item{data}{a numerical vector or matrix}
  \item{window}{an integer setting the size of the window}
  \item{step}{an integer setting the size of step between windows}
  \item{strict}{\code{TRUE} or \code{FALSE} indicating whether validation testing should be performed. Defaults to \code{FALSE}.}
}
\details{
If the input data is a vector then the function is applied to each window of the specified size, stepping through the vector by the step size. If the input data is a matrix then the function is applied to square submatrices of the specified window size, stepping through both rows and columns by the step size.}
\value{
If \code{data} is a vector, returns a numeric vector of values representing the application of the selected function within each window.

If \code{data} is a matrix, returns a numeric matrix of values representing the application of the selected function within each square window.
}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{
# vector example
x1 <- rnorm(100, sd=3)
z1 <- SlidingWindow(FUN="mean", data=x1, window=10, step=5, strict=TRUE)

# matrix example
x2 <- matrix(rnorm(10000),100,100)
z2 <- SlidingWindow(FUN="mean", data=x2, window=10, step=5, strict=TRUE)

}
\keyword{ basic stats }
\keyword{ sliding window }
