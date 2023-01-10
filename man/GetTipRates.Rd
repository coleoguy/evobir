\name{GetTipRates}
\alias{GetTipRates}
\title{Calculate the rate of evolution on the leaves of a phylogeny}
\description{This function calculates the rate of change from parent nodes to extant tips of a phylogeny.}
\usage{
GetTipRates(tree, Q, tip.states, hyper, p.mat)
}
\arguments{
  \item{tree}{an ultrametric tree of class phylo}
  \item{Q}{transition matrix containing estimated rates with column names}
  \item{tip.states}{An integer vector with a length equal to the number of species on the phylogeny. It should have values of 1 to N with N being the number of total states. Elements of the vector should match the tip names for the phylogeny.}
  \item{hyper}{logical vector of length one. TRUE indicates the model includes a binary hyperstate. Default is FALSE and indicates no binary hyperstate}
  \item{p.mat}{a probability matrix where each column represent a discrete state and each row is a species. Values in the matrix describe the probability of being in any of these states}
}

\value{A named numeric vector with the rate for each tip in the phylogeny.}
\author{Michelle M. Jonika and Heath Blackmon}
\references{http://coleoguy.github.io/}
\keyword{comparative phylogenetics, q-matrix, discrete trait}
