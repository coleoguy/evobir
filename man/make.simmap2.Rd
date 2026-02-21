\name{make.simmap2}
\alias{make.simmap2}
\title{Modified stochastic mapping which is resistant to model failure}
\description{This function is an extension of the \code{make.simmap} function of \pkg{phytools} which allows users to monitor rejections at specific edges during the simulation process and set an upper limit on the number of rejections permitted per edge. When the rejection limit is exceeded for a given edge, it is assigned a \code{"fail"} state rather than causing the entire simulation to stall. Failed edges can later be resolved using \code{\link{fix.simmap}}.}
\usage{
make.simmap2(tree, x, model="SYM", nsim=1 ,
             monitor = FALSE, rejmax = NULL, rejint = 1000000,
             ...)
}
\arguments{
  \item{tree}{a phylogenetic tree of class \code{phylo} or \code{multiPhylo}. See \code{\link[phytools]{make.simmap}} for details.}
  \item{x}{a named vector or matrix of tip states. See \code{\link[phytools]{make.simmap}} for details.}
  \item{model}{a character string or matrix describing the model. Options include \code{"ER"}, \code{"SYM"}, \code{"ARD"}, or a custom rate matrix. See \code{\link[phytools]{make.simmap}} for details.}
  \item{nsim}{number of stochastic maps to simulate. Defaults to 1.}
    \item{monitor}{logical; whether or not to print the number of rejections per edge to the console. Defaults to \code{FALSE}.}
  \item{rejmax}{integer giving the maximum number of rejections permitted before an edge is skipped. Defaults to \code{NULL} (no upper limit on rejections).}
  \item{rejint}{integer giving the interval at which the rejection count is printed to the console (when \code{monitor = TRUE}). Defaults to 1000000.}
  \item{...}{optional arguments passed to \code{\link[phytools]{make.simmap}}, such as \code{Q}, \code{pi}, \code{burnin}, \code{samplefreq}, or \code{prior}.}
}
\value{An object of class \code{simmap} (if \code{nsim = 1}) or \code{multiSimmap} (if \code{nsim > 1}). See \code{\link[phytools]{make.simmap}} for general structure.

In addition to the states present in the model, an additional state \code{"fail"} in the \code{maps} and \code{mapped.edge} elements is assigned to edges which are skipped due to exceeding the rejection limit. A \code{fail} element on the tree object lists the edge indices that failed.
}
\seealso{
  \code{\link[phytools]{make.simmap}}, \code{\link{fix.simmap}}
}
\author{Matthew Marano, Maximos Chin, and Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{
\dontrun{
library(ape)
library(phytools)
## Simulate a tree and discrete character
tree <- rtree(20)
x <- sample(c("a","b"), 20, replace = TRUE)
names(x) <- tree$tip.label

## Run stochastic mapping with rejection monitoring
result <- make.simmap2(tree, x, model = "ER", nsim = 1,
                       monitor = FALSE, rejmax = 100000)
}
}
\keyword{comparative phylogenetics}
\keyword{discrete trait}
\keyword{stochastic mapping}
