\name{GetTipRates}
\alias{GetTipRates}
\title{Calculate the rate of evolution on the leaves of a phylogeny}
\description{This function calculates the rate of change from parent nodes to extant tips of a phylogeny. It uses ancestral state reconstruction via the Mk model (implemented in \pkg{castor}) to estimate states at internal nodes, then computes the rate of change along each terminal branch as the absolute difference between the parent node state and tip state divided by the branch length.}
\usage{
GetTipRates(tree = NULL, Q = NULL, tip.states = NULL,
            hyper = FALSE, p.mat = NULL)
}
\arguments{
  \item{tree}{an ultrametric tree of class \code{phylo}}
  \item{Q}{a transition matrix containing estimated rates with column names. Column names should correspond to the actual state values.}
  \item{tip.states}{a named numeric vector with a length equal to the number of species on the phylogeny. Values should be integers from 1 to N where N is the number of total states. Names should match the tip labels of the phylogeny. Can be \code{NULL} if \code{p.mat} is provided.}
  \item{hyper}{logical; \code{TRUE} indicates the model includes a binary hyperstate (state labels contain an 'h' prefix that will be stripped). Defaults to \code{FALSE}.}
  \item{p.mat}{a probability matrix where each column represents a discrete state and each row is a species. Values describe the probability of being in each state. Row names should match tree tip labels. Used to generate \code{tip.states} when \code{tip.states} is \code{NULL}. Defaults to \code{NULL}.}
}

\value{A named numeric vector with the rate of change for each tip in the phylogeny. Names correspond to tip labels.}
\author{Michelle M. Jonika and Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{
\donttest{
## Create a simple example
library(ape)
tree <- rtree(10)
tree <- chronos(tree)
class(tree) <- "phylo"

## Define a simple 2-state Q matrix
Q <- matrix(c(-0.5, 0.5, 0.5, -0.5), 2, 2)
colnames(Q) <- rownames(Q) <- c("1", "2")

## Create tip states
tip.states <- sample(1:2, 10, replace = TRUE)
names(tip.states) <- tree$tip.label

## Calculate tip rates
rates <- GetTipRates(tree = tree, Q = Q, tip.states = tip.states, hyper = FALSE)
rates
}
}
\keyword{comparative phylogenetics}
\keyword{discrete trait}
