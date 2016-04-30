\name{StochRound}
\alias{StochRound}
\title{
Stochastic rounding}
\description{
Rounds a number up or down to a whole number with a probability equal to how close the starting value is to these whole numbers.}
\usage{
StochRound(x)
}
\arguments{
  \item{x}{a numerical value.}
}
\details{
Stochastis rounding can be helpful in simulations where we need a round number (for instance number individuals) but we calculate based on some percentage.  Using this function we will pick reasonable numbers without a bias.}
\value{
one of the two closest whole numbers is returned
}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{
StochRound(7.4)
}
\keyword{ basic stats }
\keyword{ round }
