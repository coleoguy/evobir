\name{countTrees}
\alias{countTrees}
\title{Calculate the number of times a set of topologies occur}
\description{This function counts the number of times that a set of topologies is present in a collection of trees.}
\usage{
countTrees(collection = NULL, ref = NULL, classes = T, verbose=T)
}
\arguments{
  \item{collection}{path to a collection of trees in a Newick format file}
  \item{ref}{path to a Newick format file with the topologies of interest}
  \item{classes}{if T then will return a vector with classification of each tree}
}
\details{}

\value{If classes is T returns a list with the first element being a numeric vector of the same length as the number of trees in the ref file. The elements of the returned vector correspond to the occurences of trees in the collection file that match the topologies supplied in the ref file. The second element is a vector the same length as the input tree collection and each tree is assigned a number based on the topology it matches.}
\author{Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{}
\keyword{}
