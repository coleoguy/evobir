\name{BranchTimes}
\alias{BranchTimes}
\title{
Computes branch ages}
\description{
Computes the time measured from the tip of the tree for each branch's origin.}
\usage{
BranchTimes(tree)
}
\arguments{
  \item{tree}{an ultrametric phylogeny of class phylo}
}
\details{
Returns a numeric vector of distances from the tip of the tree to the origin of each branch on the phylogeny.  This vector is in the same order as the branches are listed in the matrix tree$edge, and that the lengths are listed in tree$edge.length.}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{
data(trees.mite)
BranchTimes(tree[[1]])
}
\keyword{ phylogeny }
\keyword{ tree }
