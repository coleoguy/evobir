\name{FuzzyMatch}
\alias{FuzzyMatch}
\title{Find Close Matches in a tree and dataset}
\description{When assembling data from different sources typos can sometimes cause a loss of perfect matches between trees and datasets.  This function helps you find these close matches that can be hand curated to keep as many species as possible in your analysis.}
\usage{
FuzzyMatch(tree, data, max.dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{ a phylogenetic tree of the class "phylo".}
  \item{data}{ character vector with the names from your dataset.}
  \item{max.dist}{ This is the maximum number of characters that can differ between your tree and data and still be recognized as a close match.}
}

\value{
  A dataframe with the following rows:

  Name in data\cr
  Name in tree\cr
  Number of differences\cr
}
\references{http://coleoguy.github.io/}
\author{Heath Blackmon}
\examples{
data(hym.tree)
names <- c("Pepsis_elegans", "Plagiolepis_alluaudi", "Pheidele_lucreti",
           "Meliturgula_scriptifronsi", "Andrena_afimbriat")
FuzzyMatch(tree = hym.tree, data = names, max.dist=3)
}
\keyword{phylogenetics}
\keyword{comparative methods}
