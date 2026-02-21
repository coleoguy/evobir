\name{SampleTrees}
\alias{SampleTrees}
\title{Select a random sample of trees}
\description{This function takes as its input a large collection of trees from a program like MrBayes or BEAST, removes a user-specified proportion as burn-in, and randomly samples a desired number of trees from the post-burn-in distribution.}
\usage{SampleTrees(trees, burnin, final.number, format, prefix)}
\arguments{
  \item{trees}{path to a nexus format file containing trees that the user wants to sample from}
  \item{burnin}{the proportion of trees to remove as burn-in (e.g., 0.25 removes the first 25\%)}
  \item{final.number}{the number of trees desired in the final sample}
  \item{format}{output format: \code{"new"} for Newick format or \code{"nex"} for Nexus format}
  \item{prefix}{a text string to assign as the prefix of the new tree file name}
}
\value{
Writes the sampled trees to a file and returns an object of class \code{multiPhylo}.
}
\references{\url{http://coleoguy.github.io/}}
\author{Heath Blackmon}
\examples{
\donttest{
SampleTrees(trees = system.file("trees.nex", package = "evobiR"),
            burnin = .1, final.number = 20, format = 'new', prefix = 'sample')
}
}
\keyword{ phylogenetics }
