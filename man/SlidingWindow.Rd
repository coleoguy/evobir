\name{SlidingWindow}
\alias{SlidingWindow}
\title{
Sliding window analysis}
\description{
Applies a function within a sliding window of a numeric vector or matrix.  Both the step size and the window size can be set by the user.  For the matrix impelentation the step size and window size is constrained to be the same in both the X and Y dimensions.}
\usage{
SlidingWindow(FUN, data, window, step, strict)
}
\arguments{
  \item{FUN}{a function to be applied within each window.}
  \item{data}{a numerical vector or matrix}
  \item{window}{an integer setting the size of the window}
  \item{step}{an integer setting the size of step between windows}
  \item{strict}{TRUE or FALSE indicating whether validation testing should be performed}
}
\details{
If the input data is a vector then returns a vector of numeric values representing the application of the selected function within each window.  If the input data is a matrix then returns a matrix of numeric values representing the application fo the selected function within each window.}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{
# vector example
x1 <- rnorm(100, sd=3)
z1 <- SlidingWindow(FUN="mean", data=x2, window=10, step=5, strict=T)

# matrix example
x2 <- matrix(rnorm(10000),100,100)
z2 <- SlidingWindow(FUN="mean", data=x1, window=10, step=5, strict=T)

}
\keyword{ basic stats }
\keyword{ sliding window }
