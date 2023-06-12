\name{plot}
\alias{plot.phyloscaled}
\title{Phylogenetic visualization of heterogenity in discrete character evolution}
\description{This function provides two methods for visualizing a \code{phyloscaled} tree produced by \code{scaleTreeRates}.}
\usage{
\method{plot}{phyloscaled}(tree, method = "multiply", palette = "RdYlGn",
     edge.width = 1, cex = 1, show.tip.label = T)
}


\arguments{
  \item{tree}{a tree of class \code{phylo} and \code{phyloscaled}}
  \item{method}{a \code{string} describing the method to be used for visualization. Can either be "multiply" or "color", see \strong{Details}}
  \item{palette}{a \code{string} giving the palette to be used for edge coloring, only taken into account if \code{method = "color"}. See \strong{Details}}
  \item{edge.width}{numeric value that determines branch width}
  \item{cex}{numeric value for size}
  \item{show.tip.label}{logical indicating whether to print tip labels}
}
\details{The two plotting methods currently supported are \code{"multiply"} and \code{"color"}.

\code{multiply} takes the length of each edge on the phylogeny and multiplies it by the scalar associated with the edge before plotting the scaled tree.

\code{color} assigns a color from a diverging palette to each edge depending on it's associated scalar. Currently supported palettes are any of the \pkg{RColorBrewer} palettes or the standard viridis palette (specified by string \code{"viridis"}). Because RColorBrewer sets the maximum number of distinct colors for divergent palettes to be 11, any phylogeny which has greater than 11 unique scalar bins represented within it's edges must use a viridis palette. The ultrametric phylogeny is then plotted with each edge colored accordingly
}
\value{A plotted phylogeny with edges either multiplied or colored by their associated scalars}
\seealso{
  \code{\link{scaleTreeRates}}
}
\author{Maximos Chin and Heath Blackmon}
\references{http://coleoguy.github.io/}
\keyword{comparative phylogenetics, discrete trait, rate heterogeneity}
