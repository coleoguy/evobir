\name{ReorderData}
\alias{ReorderData}
\title{
Reorders trait data to match the order of tips in a tree}
\description{
This function takes a vector, matrix, or dataframe and reorders the data to match the order of tips in a phylo object.}
\usage{
ReorderData(tree, data, taxa.names="row names")
}
\arguments{
  \item{tree}{a phylo object}
  \item{data}{a vector, matrix, dataframe set of taxa names as present in the tree and data must match.  If data is a vector it should be a named vector.  If the data is a matrix or dataframe the taxa names may be row names or present in a column.}
  \item{taxa.names}{If taxa names are present in a column the column number should be supplied.  If taxa names are the row names the argument can be set to "row names" (default setting).  If the data is being supplied in a vector this argument is not used.}
}
\details{
Returns data in the same format as supplied but reordered to match the order of tip labels.}
\references{
\url{http://www.uta.edu/karyodb/evobiR/}}
\author{
Heath Blackmon
}

\examples{
data(finch)
phy <- finch$phy
data <- finch$data[,"wingL"]
# we create a second dataset with the order of taxa randomized
data <- sample(data, size=13)
names(data) == phy$tip.label ## order is different
data <- ReorderData(phy, data)
names(data) == phy$tip.label ## order is the same
}
\keyword{ model comparison }
\keyword{ AICc }
