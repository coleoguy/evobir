\name{fix.simmap}
\alias{fix.simmap}
\title{Fix a stochastic map with failed edges}
\description{This function facilitates automated resolution of failed edges in a modified stochastic map produced by \code{\link{make.simmap2}} through application of graph theory implemented in \pkg{igraph}.}
\usage{
fix.simmap(hists, tips, transition.matrix)
}
\arguments{
  \item{hists}{an object of class \code{simmap} or \code{multiSimmap} produced by \code{\link{make.simmap2}}}
  \item{tips}{two column dataframe with first column listing species name as shown in tree tips and second column listing model state for each species}
  \item{transition.matrix}{square matrix describing possible transitions between states}
}
\value{A object of class \code{simmap} or \code{multiSimmap}, see \code{\link{make.simmap}}. 

All edges which failed in the original stochastic maps should be resolved.
}
\seealso{
  \code{\link{make.simmap2}}
}
\author{Maximos Chin, Matthew Marano, and Heath Blackmon}
\references{http://coleoguy.github.io/}
\keyword{comparative phylogenetics, discrete trait, stochastic mapping}
