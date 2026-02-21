\name{plot}
\alias{plot.phyloscaled}
\title{Phylogenetic visualization of heterogeneity in discrete character evolution}
\description{This function provides two methods for visualizing a \code{phyloscaled} tree produced by \code{\link{scaleTreeRates}}.}
\usage{
\method{plot}{phyloscaled}(x, method = "multiply", palette = "RdYlGn",
     edge.width = 1, cex = 1, show.tip.label = T, \dots)
}


\arguments{
  \item{x}{a tree of class \code{phylo} and \code{phyloscaled}, as produced by \code{\link{scaleTreeRates}}}
  \item{method}{a string describing the method to be used for visualization. Can either be \code{"multiply"} or \code{"color"}, see \strong{Details}.}
  \item{palette}{a string giving the palette to be used for edge coloring, only taken into account if \code{method = "color"}. See \strong{Details}.}
  \item{edge.width}{numeric value that determines branch width}
  \item{cex}{numeric value for tip label size}
  \item{show.tip.label}{logical indicating whether to print tip labels}
  \item{\dots}{additional arguments passed to plotting functions}
}
\details{The two plotting methods currently supported are \code{"multiply"} and \code{"color"}.

\code{"multiply"} takes the length of each edge on the phylogeny and multiplies it by the scalar associated with the edge before plotting the scaled tree. This produces a tree where branch lengths reflect the estimated rate variation.

\code{"color"} assigns a color from a diverging palette to each edge depending on its associated scalar. Currently supported palettes are any of the \pkg{RColorBrewer} diverging palettes or the standard viridis palette (specified by string \code{"viridis"}). Because \pkg{RColorBrewer} sets the maximum number of distinct colors for divergent palettes to be 11, any phylogeny which has greater than 11 unique scalar bins represented within its edges will automatically use the viridis palette. The ultrametric phylogeny is then plotted with each edge colored accordingly.
}
\value{A plotted phylogeny with edges either multiplied or colored by their associated scalars. Called for its side effect (plotting).}
\seealso{
  \code{\link{scaleTreeRates}}
}
\author{Maximos Chin and Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{
\donttest{
## Requires a phyloscaled tree from scaleTreeRates
library(ape)
library(phytools)
tree <- rtree(20)
tip.states <- sample(1:2, 20, replace = TRUE)
names(tip.states) <- tree$tip.label
scaled.tree <- scaleTreeRates(tree, tip.states, model = "ER", nbins = 3, max.ratio = 2)
plot(scaled.tree, method = "multiply")
}
}
\keyword{comparative phylogenetics}
\keyword{discrete trait}
\keyword{rate heterogeneity}
