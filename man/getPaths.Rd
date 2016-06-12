\name{getPaths}
\alias{getPaths}
\title{
Get evolutionary pathways}
\description{
This function takes a tree of class phylo and returns all of the pathways from the root to the tips.}
\usage{
getPaths(tree, type="branch")
}
\arguments{
  \item{tree}{a phylo object}
  \item{type}{a text string either "node" or "branch"}
}
\value{
This function returns a list the elements of which are evolutionary pathways from the root to the tips of the tree.  This is returned either as a list of branches or node depending on the value provided to the argument type}

\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}
\keyword{ tree manipulation }