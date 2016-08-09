\name{isupper}
\alias{isupper}
\alias{islower}
\title{
check case}
\description{
These functions test whether a character is upper or lower case.  If applied to a string it will return True only if all characters match the test. For example both \code{isupper("Abbi")} and \code{islower("Abbi")} will return False.  These functions are vectorized and may be applied to vectors of any length.}
\usage{
isupper(x)
islowere(x)
}
\arguments{
  \item{x}{a vector of characters or strings}
}
\details{
returns a logical vector of the same length as x}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{
x <- c("A","a","Leslie", "ABBI", "beetle")
isupper(x)
islower(x)
}
\keyword{ text }
