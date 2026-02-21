\name{countTrees}
\alias{countTrees}
\title{Calculate the number of times a set of topologies occur}
\description{This function counts the number of times that a set of reference topologies is present in a collection of trees. Trees are compared using topology only (branch lengths are ignored). This is useful for summarizing posterior distributions of trees from Bayesian phylogenetic analyses.}
\usage{
countTrees(collection = NULL, ref = NULL, classes = T, verbose=T)
}
\arguments{
  \item{collection}{path to a collection of trees in a Newick format file}
  \item{ref}{path to a Newick format file with the reference topologies of interest}
  \item{classes}{logical; if \code{TRUE} then returns both topology counts and a classification vector assigning each tree to a topology. Defaults to \code{TRUE}.}
  \item{verbose}{logical; if \code{TRUE} prints progress messages. Defaults to \code{TRUE}.}
}

\value{If \code{classes} is \code{TRUE}, returns a list with two elements:
  \item{counts}{a numeric vector with length equal to the number of reference topologies. Each element gives the count of trees in the collection matching that topology.}
  \item{classification}{a numeric vector with length equal to the number of trees in the collection. Each tree is assigned the number of the reference topology it matches, or \code{NA} if it does not match any.}

If \code{classes} is \code{FALSE}, returns only the counts vector.
}
\author{Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{
\donttest{
## Create example tree files
library(ape)
## Write 50 random trees to a file
trees <- rmtree(50, 5)
write.tree(trees, file = tempfile(fileext = ".nwk"))

## Write 3 reference topologies
refs <- rmtree(3, 5)
write.tree(refs, file = tempfile(fileext = ".nwk"))
}
}
\keyword{tree topology}
