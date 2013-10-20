\name{SampleTrees}
\alias{SampleTrees}
\title{Select a random sample of trees}
\description{This function takes as its input a large collection of trees from a program like MrBayes or Beast and allows the user to select the number of randomly drawn trees they wish to retrieve}
\usage{SampleTrees(trees, burnin, final.number, format, prefix)}
\arguments{
  \item{trees}{a nexus format file containing trees that the user wants to sample from}
  \item{burnin}{the proportion of trees to remove as burnin}
  \item{final.number}{the number of trees desired}
  \item{format}{options are "new" or "nex" indicating to save the trees in newick format or nexus format}
  \item{prefix}{a text string to assing to the new treefile name
}
}
\value{
an object of the class "multiPhylo" is returned
}
\references{\url{http://www.uta.edu/karyodb/evobiR/}}
\author{Heath Blackmon}
\examples{
SampleTrees(trees = system.file("trees.nex", package = "evobiR"), 
            burnin = .1, final.number = 20, format = 'new', prefix = 'sample')}
\keyword{ phylogenetics }
