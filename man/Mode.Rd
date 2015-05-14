\name{Mode}
\alias{Mode}
\title{
Calcualtes the mode of a numeric vector}
\description{
R's base package function mode returns the type of object 'numeric', 'character' etc.  This give the option of an easy  to remember work around for that.}
\usage{
Mode(x)
}
\arguments{
  \item{x}{a numerical vector.}
}
\details{
Returns the most frequently occuring value in a vector. In the case of a tie it will return the mode which has the earliest initial occurence in the vector}
\value{
returns the most frequently occuring value in a series of numbers
}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{
Mode(c(1,2,3,4,5,6,2,5))
}
\keyword{ basic stats }
\keyword{ mode }
