\name{getDist}
\alias{getDist}
\title{
Measure distance between trait transitions in two characters}
\description{
This function takes two stochastic maps for binary traits and measures the distances between each of the 8 type of changes in two traits.}
\usage{
getDist(map1, map2, method=NULL)
}
\arguments{
  \item{map1}{a stochastic map generated with phytools}
  \item{map2}{a stochastic map generated with phytools}
}
\value{
This function returns a list the elements of which are vectors representing the distance between each type of change.}

\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}
\keyword{ tree manipulation }