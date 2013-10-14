\name{ReplTable}
\alias{ReplTable}
\alias{MaxTips}
\title{Maximize Retained Tips}
\description{These functions attempt to retain as many datapoints as possible during a comparative analysis.}
\usage{
ReplTable(tree, data.taxonomy, tree.taxonomy, levels=c("Genus", "Family"), verbose=TRUE)

MaxTips(tree, data, repl.table)
}
\arguments{
  \item{tree}{a phylogenetic tree of the class "phylo".}
  \item{data.taxonomy}{ a datframe containing taxonomic information for all species in the dataset see function \code{\link{GetTaxonomy}}}
  \item{tree.taxonomy}{ a datframe containing taxonomic information for all species in the dataset see function \code{\link{GetTaxonomy}}}
  \item{data}{ a dataframe containg the data you would like to use in your comparative analysis}
  \item{repl.table}{ the output of ReplTable and one of the inputs of MaxTips}  
  \item{levels}{currently not implemented}
  \item{verbose}{turns on reports to the screen.}

}

\details{
 The goal of these two functions is to maximize the number of tips that are retained in a comparative phylogenetic analysis.  This is done by searching for higher taxanomic level matches.  This is accomplished by keeping one example of any available genera not matched at the species level and one available member of each family not matched at the genus or species level.
 
 The taxonomy tables should have columns name of "Family", "Genus", and "Binomial". The dataframes (tree.taxonomy, data.taxonomy, and data) should all have row names correspoinding to the binomial.
 }
\value{
  ReplTable produces a dataframe that is used as one of the inputs for MaxTips.  It lists the tips on the tree that will be replaced and what taxonomic level matches have been found at.  MaxTips returns a pruned dataset and tree.
}
\references{http://www.uta.edu/karyodb/evobiR/}
\author{Heath Blackmon and Matt Pennell}
\examples{
data('mite.data')
data('mite.tree')
data('data.taxonomy')
data('tree.taxonomy')

replacement.table <- ReplTable(tree = mite.tree, data.taxonomy = data.taxonomy, tree.taxonomy = tree.taxonomy, levels=c("Genus", "Family"), verbose=TRUE)
new.data <- MaxTips(mite.tree, mite.data, replacement.table)}
\seealso{ \code{\link{GetTaxonomy}}
}
\keyword{comparative methods}
\keyword{phylogenetics}
