\name{fix.simmap}
\alias{fix.simmap}
\title{Fix a stochastic map with failed edges}
\description{This function facilitates automated resolution of failed edges in a modified stochastic map produced by \code{\link{make.simmap2}} through application of graph theory implemented in \pkg{igraph}. For each failed edge, it finds the shortest path between the start and end states using the transition matrix and assigns equal time to each transitional state along the path.}
\usage{
fix.simmap(hists, tips, transition.matrix)
}
\arguments{
  \item{hists}{an object of class \code{simmap} or \code{multiSimmap} produced by \code{\link{make.simmap2}} that contains one or more failed edges}
  \item{tips}{a two-column dataframe with the first column listing species names as shown in tree tip labels and the second column listing the known state for each species}
  \item{transition.matrix}{a square matrix describing possible transitions between states. Non-zero entries indicate that a direct transition is possible between those states.}
}
\value{An object of class \code{simmap} or \code{multiSimmap}, see \code{\link{make.simmap}}.

All edges which were assigned a \code{"fail"} state in the original stochastic maps will be resolved with inferred transitional histories.
}
\seealso{
  \code{\link{make.simmap2}}
}
\author{Maximos Chin, Matthew Marano, and Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{
\donttest{
library(ape)
library(phytools)
## Simulate a tree and discrete character
tree <- rtree(20)
x <- sample(c("a","b"), 20, replace = TRUE)
names(x) <- tree$tip.label

## Run stochastic mapping with low rejection limit
result <- make.simmap2(tree, x, model = "ER", nsim = 1,
                       rejmax = 10)

## Create tips dataframe
tips <- data.frame(species = tree$tip.label, state = x)

## Create transition matrix
trans.mat <- matrix(c(0, 1, 1, 0), 2, 2)
rownames(trans.mat) <- colnames(trans.mat) <- c("a", "b")

## Fix any failed edges
if(!is.null(result$fail)){
  fixed <- fix.simmap(result, tips, trans.mat)
}
}
}
\keyword{comparative phylogenetics}
\keyword{discrete trait}
\keyword{stochastic mapping}
