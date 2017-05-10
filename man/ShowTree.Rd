\name{ShowTree}
\alias{ShowTree}
\title{
Plots a tree and discrete tip data}
\description{
This function takes a tree of class phylo and a named vector of tip states.}
\usage{
ShowTree(tree, tip.vals, pch = 16, cols = NULL, show.tip.label = F, 
tip.cex = 1, label.cex = 1, type = "phylogram")
}
\arguments{
  \item{tree}{a phylo object}
  \item{tip.vals}{a vector of tip states in the same order as the order of tips on the tree.}
  \item{pch}{the symbol to be used for tip values}
  \item{cols}{a vector of colors, if the vector contains fewer colors than there are tip states it will be recycled.}
  \item{show.tip.label}{display taxa labels at the tip of the phylogeny}
  \item{tip.cex}{multiplier for tip marker size}
  \item{label.cex}{multiplier for taxa label size}
  \item{type}{type of tree to plot options are 'phylogram' or 'fan'}
}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}
\keyword{ comparative methods }
\keyword{ plotting phylogenies }
