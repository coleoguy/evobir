\name{getNe}
\alias{getNe}
\title{Calculate the variance effective population size}
\description{This function calculates the variance effective population size due to unequal sex ratio. Formula are available for both autosomal loci and X chromosome loci.}
\usage{
getNe(males = 50, females = 200, locus = "A")
}
\arguments{
  \item{males}{number of males}
  \item{females}{number of females}
    \item{locus}{"A" or "X" to denote the population size of interest is for autosomal locus or X chromosome locus respectively.}

}
\details{}

\value{Returns a numeric vector of length 1 which contains the variance effective population size.}
\author{Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{getNe(males=100, females=200)}
\keyword{ population size }
