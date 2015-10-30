\name{ShowTree}
\alias{ShowTree}
\title{
Plots a tree and discrete tip data}
\description{
This function takes a tree of class phylo and a named vector of tip states.}
\usage{
ShowTree(tree, tip.vals, cols=c("red", "blue"))
}
\arguments{
  \item{tree}{a phylo object}
  \item{data}{a vector of tip states in the same order as the order of tips on the tree.}
  \item{cols}{a vector of colors, if the vector contains fewer colors than there are tip states it will be recycled.}
}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}
\keyword{ comparative methods }
\keyword{ plotting phylogenies }
