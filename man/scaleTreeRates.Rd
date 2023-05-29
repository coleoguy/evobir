\name{scaleTreeRates}
\alias{scaleTreeRates}
\title{Phylogeneitc analysis  of heterogeneity in discrete character evolution}
\description{This function performs the phylogenetic methods for analysis of heterogenity in rates of discrete character evolution described in Jonika et al. (2023).}
\usage{
scaleTreeRates(tree, tip.states, model, fixedQ = NULL, 
               max.ratio = 2, nbins = 10, max.transition = 1, 
               var.start = FALSE, pi = "Fitzjohn")
}
\arguments{
  \item{tree}{a tree of class \code{phylo}}
  \item{tip.states}{a named vector of tip states for some discrete character which is associated with the phylogeny. Order can differ from order of tips on phylogeny}
  \item{model}{the model which should be used to perform likelihod calculations. This can either be a  \code{string} which can be passed to the \code{fitMk} function of \pkg{phytools} or a symmetrical transition \code{matrix} which has transitions between states categorized into some number of distinct classes}
  \item{fixedQ}{optional argument to be used when Q-matrix with pre-estimated rates is available. Deafults to \code{NULL}}
  \item{max.ratio}{\code{num} or \code{int} greater than one descirbing the maximum ratio of scalar bins to one. Defaults to 2, i.e. scalar bins range between 0.5 and 2}
  \item{nbins}{\code{int} giving the number of scalar bins above and below 1. Defaults to 10, i.e. 10 bins below 1 and 10 bins above 1 for a total of 21 bins inclusive of 1}
  \item{max.transition}{\code{int} giving the maximum number of bins which the scalar associated with an edge can differ from the scalar associated with it's parent edge. Defaults to 1}
  \item{var.start}{\code{logical} whether or not to increment scalar values at the root of the tree. If \code{TRUE}, the analysis will be iterated across all pssible root scalar values and the best tree (highest likelihood) returned. If \code{FALSE} (default), root scalar is set to one and only a single iteration is performed}
  \item{pi}{\code{string} giving method to be used for estimating prior. Takes any option which can be passed to phytool's \code{fitMk}, defaults to \code{"fitzjohn"}}
}

\value{A phylogeny of class \code{phylo} and \code{phyloscaled}. Phylogeny has all elements normally included in an object of class \code{phylo}, with an additional element:
\item{scalar}{a \code{numeric} vector of scalars equal in length to the number of edges in phylogeny. Ordering of scalars is identical to the ordering of edges}}
\seealso{
  \code{\link{plot.phylosaved}}
}
\author{Maximos Chin and Heath Blackmon}
\references{http://coleoguy.github.io/}
\keyword{comparative phylogenetics, discrete trait, rate heterogeneity}
